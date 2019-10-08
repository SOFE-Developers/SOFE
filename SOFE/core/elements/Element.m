classdef Element < SOFE
  properties
    dimension
    nV, nB
    order
    doFTuple
    conformity
    isLagrange
    quadRule
  end
  methods % constructor
    function obj = Element(dimension, nV, nB, order)
      obj.dimension = dimension;
      obj.nV = nV;
      obj.nB = nB;
      obj.order = order;
      obj.quadRule = obj.getQuadRule();
    end
    function R = getDoFEnum(obj, dim, varargin) % [orient]
      if size(obj.doFTuple,1) == 1, obj.doFTuple(2,:) = 1; end
      R = (1:prod(obj.doFTuple(:,dim+1)))';
    end
    function R = getQuadRule(obj, varargin) % [order]
      try r = varargin{:}; catch, r = max(max(2*(obj.order)),1); end
      switch obj.nV(end)
        case 1
          R{1} = GaussPoint();
        case 2
          R{2} = GaussPoint();
          R{1} = GaussInt(r);
        case 3
          R{3} = GaussPoint();
          R{2} = GaussInt(r);
          R{1} = GaussTri(r);
        case 4
          if obj.dimension == 2
            R{3} = GaussPoint();
            R{2} = GaussInt(r);
            R{1} = GaussQuad(r);
          else
            R{4} = GaussPoint();
            R{3} = GaussInt(r);
            R{2} = GaussTri(r);
            R{1} = GaussTet(r);
          end
        case 8
          R{4} = GaussPoint();
          R{3} = GaussInt(r);
          R{2} = GaussQuad(r);
          R{1} = GaussHex(r);
      end
    end
    function [Rp, Rw] = getQuadData(obj, codim)
      Rp = obj.quadRule{codim+1}.points;
      Rw = obj.quadRule{codim+1}.weights;
    end
  end
  methods(Static = true)
    function R = create(nV, dimP)
      switch nV
        case 1
          R = E0();
        case 2
          R = PpL(1,1);
        case 3
          R = PpL(2,1);
        case 4
          if dimP == 2
            R = QpL(2,1);
          else
            R = PpL(3,1);
          end
        case 8
          R = QpL(3,1);
      end
    end
  end
  methods % evaluation
    function R = evalBasis(obj, points, order) %#ok<STOUT,INUSL>
      eval(['R = obj.evalD' num2str(order) 'Basis(points);']);
    end
  end
  methods % get
    function R = isSimplex(obj)
      nD = obj.dimension;
      R = (obj.nV(nD)==3 && nD==2) || (obj.nV(nD)==4 && nD==3);
    end
    function R = getNC(obj)
      R = size(obj.evalBasis(zeros(1,obj.dimension),0),3);
    end
    function R = getNEntSub(obj, dim)
      switch dim
        case 0
          R = 1;
        case 1
          R = [2 1];
        case 2
          if obj.isSimplex
            R = [3 3 1];
          else
            R = [4 4 1];
          end
        case 3
          if obj.isSimplex
            R = [4 6 4 1];
          else
            R = [8 12 6 1];
          end
      end
    end
  end
  methods % display
    function show(obj, dim, i)
      if i > obj.nB(dim)
        error('XXX: There are only %d basis functions in dim %d!\n', obj.nB(dim), dim);
      end
      gridFine = linspace(0,1,300)';
%       gridCoarse = linspace(0,1,20)';
      clf
      switch dim
        case 1
          Z = obj.evalBasis(gridFine, 0);
          plot(gridFine,Z(i,:));
        case 2
          [Xf,Yf] = meshgrid(gridFine, gridFine);
%           [Xc,Yc] = meshgrid(gridCoarse, gridCoarse);
          Z = obj.evalBasis([Xf(:) Yf(:)], 0); % nBxnPxnC
          if obj.isSimplex()
            Z(:,Xf+Yf>1) = NaN;
          end
          surf(Xf,Yf,reshape(Z(i,:),size(Xf))); 
          shading interp; view(0,90)
        case 3
          N = 30;
          [X,Y,Z] = meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));
          if obj.isSimplex()
            I = X+Y+Z<=1;
            X = X(I); Y = Y(I); Z = Z(I);
          end
          P = [X(:) Y(:) Z(:)];
          W = obj.evalBasis(P, 0);
          plot3k(P, 'ColorData', W(i,:)');
      end
    end
    function quiver2D(obj, i, varargin)
      if i > obj.nB(2)
        error('XXX: There are only %d basis functions in dim %d!\n', obj.nB(2), 2);
      end
      LS = linspace(0,1,20)';
      [X,Y] = meshgrid(LS, LS);
      Z = obj.evalBasis([X(:) Y(:)], 0); % nBxnPxnC
      if obj.nV(2) == 3
        Z(:,X+Y>1) = NaN;
      end
      quiver(X,Y,reshape(Z(i,:,1),size(X)),reshape(Z(i,:,2),size(X)), 1, varargin{:});
      axis([0 1 0 1]); axis equal
    end
    function quiver3D(obj, i, varargin)
      if i > obj.nB(3)
        error('XXX: There are only %d basis functions in dim %d!\n', obj.nB(3), 3);
      end
      LS = linspace(0,1,10)';
      [X,Y,Z] = meshgrid(LS, LS, LS);
      W = obj.evalBasis([X(:) Y(:) Z(:)], 0); % nBxnPxnC
      if obj.nV(2) == 3
        W(:,X+Y+Z>1) = NaN;
        plot3([0 1 0 0 0 1 0 0],[0 0 1 0 0 0 0 1],[0 0 0 0 1 0 1 0]);
        hold on
      end
      quiver3(X,Y,Z,reshape(W(i,:,1),size(X)), ...
                   reshape(W(i,:,2),size(X)), ...
                   reshape(W(i,:,3),size(X)), varargin{:});
      hold off
      axis([0 1 0 1 0 1]); axis equal
    end
  end
  methods % prolongation
    function R = getProlongationCoeff(obj)
      [P0,D0] = obj.getLagrangeFunctionals(0);
      F0 = sum(bsxfun(@times,obj.evalBasis(P0,0),permute(D0,[3 1 2])),3)';
      R = cell(1,1,2^obj.dimension);
      for k = 1:numel(R)
        [P1,D1] = obj.getLagrangeFunctionals(k);
        R{k} = F0 \ sum(bsxfun(@times,obj.evalBasis(P1,0),permute(D1,[3 1 2])),3)';
        R{k}(abs(R{k})<1e-12) = 0;
      end
      R = cell2mat(R);
    end
  end
  methods % Lagrange functionals
    function [R,D] = getLagrangeFunctionals(obj, childNr, varargin)
      [R,D] = obj.getLagrangePoints(obj.dimension, obj.order, varargin{:});
      if childNr>0
        if obj.isSimplex()
          switch obj.dimension
            case 1
              assert(0, 'No simplex!');
            case 2
              m = Mesh([0 0; 1 0; 0 1], [1 2 3]);
            case 3
              m = Mesh([0 0 0; 1 0 0; 0 1 0; 0 0 1], [1 2 3 4]);
          end
        else
          switch obj.dimension
            case 1
              m = Mesh([0;1], [1 2]);
            case 2
              m = Mesh([0 0; 1 0; 0 1; 1 1], [1 2 3 4]);
            case 3
              m = Mesh([0 0 0; 1 0 0; 0 1 0; 1 1 0; ...
                        0 0 1; 1 0 1; 0 1 1; 1 1 1], [1 2 3 4 5 6 7 8]);
          end
        end
        m.nodes = m.topology.uniformRefine()*m.nodes;
        R = permute(m.evalReferenceMap(R, 0, childNr),[2 3 1]); % nPxnD
        switch obj.conformity
          case 'H1'
            % do nothing
          case 'HRot'
            DPhi = permute(m.evalReferenceMap(R,1,childNr), [2 3 4 1]); % nPxnwxnD
            D = sum(bsxfun(@times, DPhi, permute(D,[1 3 2])), 3); % nPxnw
          case 'HDiv'
            [~, DPhiInv, J] = m.evalTrafoInfo(R,childNr); % nExnPxnDxnW, nExnP
            DPhiInv = permute(DPhiInv, [2 4 3 1]); % nPxnWxnD ! transposed
            J = permute(J, [2 1]); % nP
            D = sum(bsxfun(@times, DPhiInv, permute(D,[1 3 2])), 3); % nPxnw
            D = bsxfun(@times, D, J); % nPxnw
        end
      end
    end
    function [R,D] = getLagrangeFunctionals_(obj, childNr)
      [R,D] = obj.getLagrangePoints(obj.dimension, obj.order);
      switch obj.dimension
        case 1
          error('TODO');
        case 2
          if obj.isSimplex
            switch childNr
              case 0
              case 1
                R = 0.5*R;
              case 2
                R = 0.5*R;
                R(:,1) = R(:,1) + 0.5;
              case 3
                R = 0.5*R;
                R(:,2) = R(:,2) + 0.5;
              case 4
                R = -0.5*R;
                R = R + 0.5;
                if ~strcmp(obj.conformity, 'H1')
                  D = -D;
                end
            end
          else
            switch childNr
              case 0
              case 1
                R = 0.5*R;
              case 2
                R = 0.5*R;
                R(:,1) = R(:,1) + 0.5;
              case 3
                R = 0.5*R;
                R(:,2) = R(:,2) + 0.5;
              case 4
                R = 0.5*R;
                R = R + 0.5;
            end
          end
        case 3
          if obj.isSimplex
            error('TODO');
          else
            switch childNr
              case 0
              case 1
                R = 0.5*R;
              case 2
                R = 0.5*R;
                R(:,1) = R(:,1) + 0.5;
              case 3
                R = 0.5*R;
                R(:,2) = R(:,2) + 0.5;
              case 4
                R = 0.5*R;
                R(:,[1 2]) = R(:,[1 2]) + 0.5;
              case 5
                R = 0.5*R;
                R(:,3) = R(:,3) + 0.5;
              case 6
                R = 0.5*R;
                R(:,[1 3]) = R(:,[1 3]) + 0.5;
              case 7
                R = 0.5*R;
                R(:,[2 3]) = R(:,[2 3]) + 0.5;
              case 8
                R = 0.5*R;
                R = R + 0.5;
            end
          end
      end
      if childNr>0 && ~strcmp(obj.conformity, 'H1')
        if obj.dimension==3 && strcmp(obj.conformity, 'HDiv') 
          D = 0.25*D;
        else
          D = 0.5*D;
        end
      end
    end
  end
  methods % for matrix free plugin
    function R = getInnerDoFKey(obj, dim)
      switch dim
        case 0
          R = obj.getInnerDoFKeyNode();
        case 1
          R = obj.getInnerDoFKeyEdge();
        case 2
          R = obj.getInnerDoFKeyFace();
      end
      if isempty(R)
        error("Inner DoFKey must not be empty!");
      end
    end
    function R = getInnerDoFKeyH(obj, dim)
      switch dim
        case 0
          R = obj.getInnerDoFKeyNodeH();
        case 1
          R = obj.getInnerDoFKeyEdgeH();
        case 2
          R = obj.getInnerDoFKeyFaceH();
        case 3
          R = obj.getInnerDoFKeyCellH();
      end
      if isempty(R)
        error("Inner DoFKey must not be empty!");
      end
    end
    function R = getInnerDoFKeyNode(obj)
      dTp = obj.doFTuple;
      if size(dTp,1)>1
        dTp = prod(obj.doFTuple);
      end
      switch obj.dimension
        case 1
          assert(false, 'TODO');
        case 2
          R = kron(eye(obj.nV(2)), ones(1,dTp(1)));
          if obj.isSimplex()
            S = [1 0 1; 1 1 0; 0 1 1];
          else
            S = [1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1];
          end
          R = [R, kron(S, ones(1,dTp(2)))];
          R = [R, ones(size(R,1), dTp(3))];
        case 3
          R = kron(eye(obj.nV(3)), ones(1,dTp(1)));
          if obj.isSimplex()
            S = zeros(4, 6);
            S(1,[1 3 4]) = 1; S(2,[1 2 5]) = 1; S(3,[2 3 6]) = 1; S(4,[4 5 6]) = 1;
            T = zeros(4, 4);
            T(1,[1 2 4]) = 1; T(2,[1 2 3]) = 1; T(3,[1 3 4]) = 1; T(4,[2 3 4]) = 1;
          else
            S = zeros(8, 12);
            S(1,[1 5 9]) = 1; S(2,[1 7 10]) = 1; S(3,[2 5 11]) = 1; S(4,[2 7 12]) = 1;
            S(5,[3 6 9]) = 1; S(6,[3 8 10]) = 1; S(7,[4 6 11]) = 1; S(8,[4 8 12]) = 1;
            T = zeros(8, 6);
            T(1,[1 3 5]) = 1; T(2,[1 4 5]) = 1; T(3,[1 3 6]) = 1; T(4,[1 4 6]) = 1;
            T(5,[2 3 5]) = 1; T(6,[2 4 5]) = 1; T(7,[2 3 6]) = 1; T(8,[2 4 6]) = 1;
          end
          R = [R, kron(S, ones(1,dTp(2)))];
          R = [R, kron(T, ones(1,dTp(3)))];
          R = [R, ones(size(R,1), dTp(4))];
      end
      R = R'; % nBxnType
    end
    function R = getInnerDoFKeyEdge(obj)
      dTp = obj.doFTuple;
      if size(dTp,1)>1
        dTp = prod(obj.doFTuple);
      end
      switch obj.dimension
        case 2
          R = kron(eye(obj.nV(2)), ones(1,dTp(2)));
          R = [R, ones(size(R,1), dTp(3))];
        case 3
          if obj.isSimplex()
            R = kron(eye(6), ones(1,dTp(2)));
            S = zeros(6,4);
            S(1,[1 2]) = 1; S(2,[1 3]) = 1; S(3,[1 4]) = 1;
            S(4,[2 4]) = 1; S(5,[2 3]) = 1; S(6,[3 4]) = 1;
          else
            R = kron(eye(12), ones(1,dTp(2)));
            S = zeros(12,6);
            S(1,[1 5]) = 1; S(2,[1 6]) = 1; S(3,[2 5]) = 1; S(4,[2 6]) = 1;
            S(5,[1 3]) = 1; S(6,[2 3]) = 1; S(7,[1 4]) = 1; S(8,[2 4]) = 1;
            S(9,[3 5]) = 1; S(10,[4 5]) = 1; S(11,[3 6]) = 1; S(12,[4 6]) = 1;
          end
          R = [R, kron(S, ones(1,dTp(3)))];
          R = [R, ones(size(R,1), dTp(4))];
      end
      R = R'; % nBxnType
    end
    function R = getInnerDoFKeyFace(obj)
      assert(obj.dimension==3);
      dTp = obj.doFTuple;
      if size(dTp,1)>1
        dTp = prod(obj.doFTuple);
      end
      if obj.isSimplex()
        R = kron(eye(4), ones(1,dTp(3)));
        R = [R, ones(size(R,1), dTp(4))];
      else
        R = kron(eye(6), ones(1,dTp(3)));
        R = [R, ones(size(R,1), dTp(4))];
      end
      R = R'; % nBxnType
    end
    function R = getInnerDoFKeyNodeH(obj)
      dTp = obj.doFTuple;
      if size(dTp,1)>1
        dTp = prod(obj.doFTuple);
      end
      nESub = obj.getNEntSub(obj.dimension);
      switch obj.dimension
        case 1
          assert(false, 'TODO');
        case 2
          R = kron(eye(nESub(1)), ones(1,dTp(1)));
          R = [R, zeros(size(R,1), nESub(2)*dTp(2)+dTp(3))];
        case 3
          R = kron(eye(nESub(1)), ones(1,dTp(1)));
          R = [R, zeros(size(R,1), nESub(2)*dTp(2)+nESub(3)*dTp(3)+nESub(4)*dTp(4))];
      end
      R = R'; % nBxnType
    end
    function R = getInnerDoFKeyEdgeH(obj)
      dTp = obj.doFTuple;
      if size(dTp,1)>1
        dTp = prod(obj.doFTuple);
      end
      nESub = obj.getNEntSub(obj.dimension);
      switch obj.dimension
        case 2
          R = zeros(nESub(2), nESub(1)*dTp(1));
          R = [R, kron(eye(nESub(2)), ones(1,dTp(2)))];
          R = [R, zeros(size(R,1), nESub(3)*dTp(3))];
        case 3
          R = zeros(nESub(2), nESub(1)*dTp(1));
          R = [R, kron(eye(nESub(2)), ones(1,dTp(2)))];
          R = [R, zeros(size(R,1), nESub(3)*dTp(3)+nESub(4)*dTp(4))];
      end
      R = R'; % nBxnType
    end
    function R = getInnerDoFKeyFaceH(obj)
      dTp = obj.doFTuple;
      if size(dTp,1)>1
        dTp = prod(obj.doFTuple);
      end
      nESub = obj.getNEntSub(obj.dimension);
      switch obj.dimension
        case 2
          R = zeros(nESub(3), nESub(1)*dTp(1)+nESub(2)*dTp(2));
          R = [R, ones(size(R,1), nESub(3)*dTp(3))];
        case 3
          R = zeros(nESub(3), nESub(1)*dTp(1)+nESub(2)*dTp(2));
          R = [R, kron(eye(nESub(3)), ones(1,dTp(3)))];
          R = [R, zeros(size(R,1), nESub(4)*dTp(4))];
      end
      R = R'; % nBxnType
    end
    function R = getInnerDoFKeyCellH(obj)
      dTp = obj.doFTuple;
      if size(dTp,1)>1
        dTp = prod(obj.doFTuple);
      end
      nESub = obj.getNEntSub(obj.dimension);
      switch obj.dimension
        case 2
          assert(0, 'Not allowed!');
        case 3
          R = zeros(nESub(4), nESub(1)*dTp(1)+nESub(2)*dTp(2)+nESub(3)*dTp(3));
          R = [R, ones(size(R,1), nESub(4)*dTp(4))];
      end
      R = R'; % nBxnType
    end
  end
  methods(Static = true)
    function R = getLegendreCoefficients(N, varargin) % [order]
      if ~isempty(varargin), order = varargin{1}; else, order = 0; end
      if N==0
        R = 1;
      elseif N==1
        R = [0 1;1 0];
      else
        R = zeros(N+1);
        R(1,end) = 1;
        R(2,end-1) = 1;
        for n = 3:N+1
          R(n,:) = (-(n-2)*R(n-2,:) + (2*n-3)*[R(n-1,2:end) 0])/(n-1);
        end
      end
      if order==-1
        R = bsxfun(@rdivide,R,(N+1:-1:1)); % R = R./(N+1:-1:1);
        R = [R -R*(-1).^(size(R,2):-1:1)'];
      else
        for k = 1:order
          R = bsxfun(@times,R(:,1:end-1),(N+1-k:-1:1));
        end
      end
    end
    function R = evalLegendre(x, N, order)
      R = zeros(numel(x),numel(N));
      for i = 1:numel(N)
        n = N(i);
        if n < order
          R(:,i) = zeros(size(x));
        else
           c = Element.getLegendreCoefficients(n, order);
           R(:,i) = polyval(c(end,:), x);
        end
      end
    end
    function R = evalShapeFunctions(x, N, k)
      R = zeros(numel(x),numel(N));
      for i = 1:numel(N)
        n = N(i);
        if n==0
          switch k
            case 0
              R(:,i) = 0.5*(1-x);
            case 1
              R(:,i) = -0.5.*ones(size(x));
          end
        end
        if n==1
          switch k
            case 0
              R(:,i) = 0.5*(1+x);
            case 1
              R(:,i) = 0.5.*ones(size(x));
          end
        end
        if n>1 && k<n+1
          c = sqrt((2*n-1)/2)*Element.getLegendreCoefficients(n-1, k-1);
          R(:,i) = polyval(c(end,:), x);
        end
      end
    end
    function R = evalReducedLegendre(x, N, k)
      R = zeros(numel(x), N-1);
      for n = N:-1:2
        c = 2*cumsum(sqrt((2*n-1)/2)*Element.getLegendreCoefficients(n-1,-1),2);
        c = 2*cumsum(c(:,1:end-1).*(-1).^(2:n+1),2);
        c = c(:,1:end-1).*(-1).^(1:n-1);
        %
        switch k
          case 0
            c = c(end,1:end);
          case 1
            c = c(end,1:end-1).*(n-2:-1:1);
          case 2
            c = (c(end,1:end-2).*(n-2:-1:2)).*(n-3:-1:1);
          case 3
            c = ((c(end,1:end-3).*(n-2:-1:3)).*(n-3:-1:2)).*(n-4:-1:1);
        end
        R(:,n-1) = polyval(c, x);
      end
    end
  end
end

