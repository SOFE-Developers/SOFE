classdef Mesh < SOFE
  properties
    element
    topology
  end
  methods % constructor & more
    function obj = Mesh(nodes, elem, varargin) % [dimP]
      if ~isempty(elem)
        obj.initMesh(nodes, elem, varargin{:});
      end
    end
    function initMesh(obj, nodes, elem, varargin)
      if nargin > 3, dimP = varargin{:}; else, dimP = size(nodes, 2); end
      obj.element = obj.getShapeElement(size(elem,2), dimP);
      obj.topology = obj.getTopology(nodes, elem, dimP);
    end
  end
  methods % reference map
    function R = evalReferenceMap(obj, points, order, varargin) % [I]
      I = ':'; if nargin > 3, I = varargin{1}; end
      if isempty(points)
        R = permute(obj.topology.nodes(I,:), [1 3 2]); % nEx1xnW
        return;
      end
      if iscell(points)
        I = points{2}; points = points{1};
        pVecN = [1 3 4 5 2]; pVecB = [2 3 4 5 1];
      else
        pVecN = [1 4 3 5 2]; pVecB = [5 2 3 4 1];
      end
      entity = obj.topology.getEntity(size(points,2), I); % nExnB
      B = obj.element.evalBasis(points, order); % nBxnPx1[xnD]
      N = reshape(obj.topology.nodes(entity(:),:),[],size(B,1),obj.topology.dimW); % nExnBxnW
      R = sum(bsxfun(@times, permute(N,pVecN), permute(B,pVecB)),5); % nExnPxnW[xnD] or nExnW[xnD]
    end
    function [R, invR, jacR] = evalTrafoInfo(obj, points, varargin) % [I]
      R = obj.evalReferenceMap(points, 1, varargin{:}); % nExnPxnWxnD or nExnWxnD
      nW = obj.topology.dimW;
      if iscell(points)
        nD = size(points{1}, 2); nP = -1;
      else
        [nP,nD] = size(points);
        R = reshape(R,[],nW,nD); % (nE*nP)xnWxnD
      end
      switch nD
        case 1 % line
          switch nW
            case 1
              invR = 1.0;
              jacR = R; 
              invR = bsxfun(@rdivide, invR, jacR);
            case {2,3}
              RTR = sum(bsxfun(@times, permute(R,[1 3 4 2]), permute(R,[1 4 3 2])),4); % (nE*nP)x1x1
              invR = 1;
              jacR = RTR; % (nE*nP)
              invR = bsxfun(@rdivide, invR, jacR);
              invR = sum(bsxfun(@times, permute(invR,[1 2 4 3]), permute(R,[1 4 2 3])),4); % (nE*nP)x1xnW
              jacR = sqrt(jacR); % (nE*nP)
          end
        case 2 % surface
          switch nW
            case 2
              invR = -R;
              invR(:,1,1) = R(:,2,2);
              invR(:,2,2) = R(:,1,1); % (nE*nP)x2x2
              jacR = R(:,1,1).*R(:,2,2) - R(:,1,2).*R(:,2,1); 
              invR = bsxfun(@rdivide, invR, jacR);
            case 3 % pseudo inverse
              RTR = sum(bsxfun(@times, permute(R,[1 3 4 2]), permute(R,[1 4 3 2])),4); % (nE*nP)x3x2
              invR = -RTR;
              invR(:,1,1) = RTR(:,2,2);
              invR(:,2,2) = RTR(:,1,1); % (nE*nP)x2x3
              jacR = RTR(:,1,1).*RTR(:,2,2) - RTR(:,1,2).*RTR(:,2,1); % (nE*nP)
              invR = bsxfun(@rdivide, invR, jacR);
              invR = sum(bsxfun(@times, permute(invR,[1 2 4 3]), permute(R,[1 4 2 3])),4); % (nE*nP)x2x3
              jacR = sqrt(jacR); % (nE*nP)
          end
        case 3 % volume
          invR(:,1,1) =   R(:,2,2).*R(:,3,3) - R(:,2,3).*R(:,3,2);
          invR(:,2,1) = -(R(:,2,1).*R(:,3,3) - R(:,3,1).*R(:,2,3));
          invR(:,3,1) =   R(:,2,1).*R(:,3,2) - R(:,2,2).*R(:,3,1);
          invR(:,1,2) = -(R(:,1,2).*R(:,3,3) - R(:,1,3).*R(:,3,2));
          invR(:,2,2) =   R(:,1,1).*R(:,3,3) - R(:,1,3).*R(:,3,1);
          invR(:,3,2) = -(R(:,1,1).*R(:,3,2) - R(:,1,2).*R(:,3,1));
          invR(:,1,3) =   R(:,1,2).*R(:,2,3) - R(:,1,3).*R(:,2,2);
          invR(:,2,3) = -(R(:,1,1).*R(:,2,3) - R(:,1,3).*R(:,2,1));
          invR(:,3,3) =   R(:,1,1).*R(:,2,2) - R(:,1,2).*R(:,2,1); % (nE*nP)x3x3
          jacR = (R(:,1,1).*R(:,2,2).*R(:,3,3) + ...
                  R(:,1,2).*R(:,2,3).*R(:,3,1) + ...
                  R(:,1,3).*R(:,2,1).*R(:,3,2) - ...
                  R(:,1,1).*R(:,2,3).*R(:,3,2) - ...
                  R(:,1,2).*R(:,2,1).*R(:,3,3) - ...
                  R(:,1,3).*R(:,2,2).*R(:,3,1)); % (nE*nP)
          invR = bsxfun(@rdivide, invR, jacR);
      end
      if nP > 0
        R = reshape(R,[],nP,nW,nD); % nExnPxnDxnW
        invR = reshape(invR,[],nP,nD,nW); % nExnPxnDxnW
        jacR = reshape(jacR,[],nP); % nExnP
      end
    end
    function R = evalInversReferenceMap(obj, points, varargin) % [output flag]
      armijo = false; armijoMax = 5; notFMax = 2; out = false;
      C = obj.topology.globalSearcher.findCandidates(points);
      [nP, nC] = size(C);
      H = zeros(nP,1); L = zeros(size(points)); Ic = (1:nP)';
      for i = 1:nC
        if isempty(Ic), break; end
        Ic(C(Ic,i)==0) = [];
        In = (1:numel(Ic))';
        InotF = zeros(size(In));
        pLoc = repmat(obj.topology.getCenterLoc(), numel(Ic),1);
        for n = 1:10
          pLocN = pLoc(In,:);
          Phi = obj.evalReferenceMap({pLocN, C(Ic(In),i)}, 0);
          F = points(Ic(In),:) - Phi;
          % test for convergence
          InotF = InotF + (~obj.topology.isFeasible(pLocN));
          normF = sum(F.^2,2); del = normF<1e-12 | InotF>notFMax;
          In(del) = []; F(del,:) = [];
          normF(del) = []; InotF(del) = []; pLocN(del,:) = [];
          if isempty(In), break; end
          % Newton step
          [~, DPhiInv] = obj.evalTrafoInfo({pLocN, C(Ic(In),i)});
          delta = sum(bsxfun(@times, DPhiInv, permute(F,[1 3 2])), 3);
          if ~armijo
            pLoc(In,:) = pLocN + delta;
          else % Armijo control
            step = ones(size(In));
            Is = (1:numel(step))';
            pLocTmp = pLocN;
            cnt = 0;
            while ~isempty(Is) && cnt < armijoMax
              cnt = cnt+1;
              pLocTmp(Is,:) = pLocN(Is,:) + bsxfun(@times, step(Is), delta(Is,:));
              PhiTmp = obj.evalReferenceMap({pLocTmp(Is,:), C(Ic(In(Is)),i)}, 0);
              normFTmp = sum((points(Ic(In(Is),:),:) - PhiTmp).^2,2);
              Is(normFTmp <= normF(Is)) = [];
              step = step/2;
            end
            pLoc(In,:) = pLocTmp;
          end
          if out, fprintf('Cand=%d(#points:%d), nNewton=%d\n',i,numel(Ic),n); end
        end
        I = obj.topology.isFeasible(pLoc);
        H(Ic(I)) = C(Ic(I),i); L(Ic(I),:) = pLoc(I,:);
        Ic = Ic(~I);
      end
      R = {L,H};
    end
  end
  methods % evaluation
    function R = evalFunction(obj, F, points, U, D, varargin) % [I]
      I = ':'; if nargin > 5, I = varargin{1}; end
      P = obj.evalReferenceMap(points, 0, I); % nExnPxnW
      [nE, nP, nD] = size(P);
      P = reshape(P, nE*nP, nD); % (nE*nP)xnW
      switch nargin(F)
        case 1 % F(x)
          R = reshape(F(P), nE, nP, []);
        case 2 % F(x,U)
          if ~iscell(U); U = {U}; end
          for i = 1:numel(U)
            U{i} = reshape(U{i}(I,:), [], size(U{i},3));
          end
          R = reshape(F(P, U), nE, nP, []); % nExnPxnC
        case 3 % F(x,U,D)
          if ~iscell(U); U = {U}; end % nExnPxnC
          if ~iscell(D); D = {D}; end % nExnPxnCxnD
          sz = size(D{1});
          for i = 1:numel(U)
            U{i} = reshape(U{i}(I,:), [], sz(3));
            D{i} = reshape(D{i}(I,:), [], sz(3), sz(4));
          end
          R = reshape(F(P, U, D), nE, nP, []); % nExnPxnCxnD
      end
    end
    function [R, RVec] = integrate(obj, func, quadRule, varargin)
      if ~isreal(func)
        func = obj.evalFunction(func, quadRule.points, [], [], varargin{:});
      end
      [~,~,trafo] = obj.evalTrafoInfo(quadRule.points,varargin{:});
      RVec = (func.*abs(trafo))*quadRule.weights;
      R =  sum(RVec);
    end
  end
  methods % refinement
    function uniformRefine(obj, N)
      for i = 1:N
        obj.topology.uniformRefine();
      end
      obj.topology.notifyObservers();
    end
  end
  methods % display
    function show(obj, type, varargin) % [type]
      obj.topology.show(varargin{:}); axis equal
      if nargin < 2 || isempty(type)
        return
      end
      hold on;
      for i = 1:numel(type)
        obj.topology.showEntity(str2double(type(i)));
      end
      hold off
      axis equal
    end
  end
  methods(Static = true)
    function R = getTopology(nodes, elem, dimP)
      switch size(elem, 2)
        case 2
          R = MeshTopologyInt(nodes, elem, dimP);
        case 3
          R = MeshTopologyTri(nodes, elem, dimP);
        case 4
          if dimP == 2
            R = MeshTopologyQuad(nodes, elem, dimP);
          else
            R = MeshTopologyTet(nodes, elem, dimP);
          end
        case 8
          R = MeshTopologyHex(nodes, elem, dimP);
      end
    end
    function R = getShapeElement(N, dimP)
      switch N
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
    function [nodes, elem] = getTensorProductMesh(grid, varargin) % [isTri]
      switch numel(grid)
        case 1
          nodes = grid{1};
          if size(nodes,1) == 1
            nodes = nodes';
          end
          m = numel(nodes);
          elem = [1:(m-1); 2:m]';
        case 2
          m = numel(grid{1});
          n = numel(grid{2});
          nodes = [kron(ones(1,n),grid{1}); ...
                   kron(grid{2}, ones(1,m))]';
          elem = [1:m*(n-1)-1; 2:m*(n-1)];
          elem = [elem; elem+m]';
          elem(m:m:end,:) = [];
          if nargin > 1
            switch varargin{1}
              case 1
                elem = [elem(:,[1 2 3]); elem(:,[4 3 2])];
              case 2
                mid = permute(sum(reshape(nodes(elem(:),:),[],4,2),2)/4,[1 3 2]);
                elem = [elem size(nodes,1)+(1:size(elem,1))'];
                nodes = [nodes; mid];
                elem = [elem(:,[5 1 2]); elem(:,[5 2 4]); ...
                        elem(:,[5 4 3]); elem(:,[5 3 1])];
            end
          end
        case 3
          m = numel(grid{1});
          n = numel(grid{2});
          p = numel(grid{3});
          nodes = [kron(ones(1,p),kron(ones(1,n),grid{1}))', ...
                   kron(ones(1,p),kron(grid{2}, ones(1,m)))', ...
                   kron(grid{3},kron(ones(1,n), ones(1,m)))'];
          elem = [1:m*n*(p-1)-m-1; 2:m*n*(p-1)-m]';
          elem = [elem elem+m];
          elem = [elem, elem+m*n];
          delete = bsxfun(@(x,y)x+y,m*(n-1)+1:m*n-1,((0:(p-3))*m*n)');
          elem([m:m:m*n*(p-1)-m-1 delete(:)'],:) = [];
          if nargin > 1 && varargin{1}
            elem = [elem(:,[5 6 1 8]); elem(:,[5 1 7 8]); elem(:,[2 1 6 8]); ...
                    elem(:,[3 7 1 8]); elem(:,[4 2 8 1]); elem(:,[3 1 4 8])];
          end
      end
    end
  end
end