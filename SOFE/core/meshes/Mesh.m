classdef Mesh < SOFE
  properties
    element
    topology
    nodes
    dimW
    globalSearcher
  end
  methods % constructor & globalsearcher
    function obj = Mesh(nodes, elem, varargin) % [dimP]
      if ~isempty(elem)
        obj.dimW = size(nodes,2);
        obj.nodes = nodes;
        if nargin > 2, dimP = varargin{1}; else, dimP = size(nodes, 2); end
        obj.element = obj.getShapeElement(size(elem,2), dimP);
        obj.topology = obj.getTopology(nodes, elem, dimP);
      end
    end
    function R = getGlobalSearcher(obj)
      if isempty(obj.globalSearcher)
        obj.globalSearcher = GlobalSearcher(obj);
      end
      R = obj.globalSearcher;
    end
  end
  methods % obj is observed
    function notify(obj)
      obj.globalSearcher = [];
      obj.notifyObservers();
    end
  end
  methods % reference map
    function R = evalReferenceMap(obj, points, order, varargin) % [I]
      I = ':'; if nargin > 3, I = varargin{1}; end
      if isempty(points)
        R = permute(obj.nodes(I,:), [1 3 2]); % nEx1xnW
        return;
      end
      if iscell(points)
        I = points{2}; points = points{1};
        if isempty(points), R = []; return; end
        pVecN = [1 3 4 5 2]; pVecB = [2 3 4 5 1];
      else
        pVecN = [1 4 3 5 2]; pVecB = [5 2 3 4 1];
      end
      B = obj.element.evalBasis(points, order); % nBxnPx1[xnD]
      entity = obj.topology.getEntity(size(points,2), I); % nExnB
      N = reshape(obj.nodes(entity(:),:),[],size(B,1),size(obj.nodes,2)); % nExnBxnW
      R = sum(bsxfun(@times, permute(N,pVecN), permute(B,pVecB)),5); % nExnPxnW[xnD] or nExnW[xnD]
    end
    function [R, invR, jacR] = evalTrafoInfo(obj, points, varargin) % [I]
      R = obj.evalReferenceMap(points, 1, varargin{:}); % nExnPxnWxnD or nExnWxnD
      nW = obj.dimW;
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
      out = false;
      armijo = false; armijoMax = 100;
      notFblMax = 3; 
      TOLF = 1e-14;
      TOLREF = 1e-12;
      %
      gs = obj.getGlobalSearcher();
      C = gs.findCandidates(points);
      [nP, nC] = size(C);
      H = zeros(nP,1); L = zeros(size(points)); Ic = (1:nP)';
      for i = 1:nC
        Ic(C(Ic,i)==0) = [];
        if isempty(Ic), break; end
        In = (1:numel(Ic))';
        cntNotFbl = zeros(size(In));
        pLoc = repmat(obj.topology.getCenterLoc(), numel(Ic),1);
        maxIt = 10;
        for n = 1:maxIt
          pLocN = pLoc(In,:);
          % test for convergence
          F = points(Ic(In),:) - obj.evalReferenceMap({pLocN, C(Ic(In),i)},0);
          cntNotFbl = cntNotFbl + ~obj.topology.isFeasible(pLocN, TOLREF);
          normF = sum(F.^2,2);
          del = normF<TOLF^2 | cntNotFbl>notFblMax;
          In(del) = []; cntNotFbl(del) = []; pLocN(del,:) = [];
          F(del,:) = []; normF(del) = [];
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
            if cnt==armijoMax, warning('!armijoMax reached!'); end
            pLoc(In,:) = pLocTmp;
          end
          if out, fprintf('Cand=%d(#points:%d), nNewton=%d\n',i,numel(Ic),n); end %#ok<UNRCH>
          if n==maxIt
            pLoc(In,:) = Inf;
          end
        end
        I = obj.topology.isFeasible(pLoc, TOLREF);
        H(Ic(I)) = C(Ic(I),i); L(Ic(I),:) = pLoc(I,:);
        Ic = Ic(~I);
      end
      R = {L,H};
    end
  end
  methods % evaluation
    function R = evalFunction(obj, F, points, S, varargin) % [I]
      I = ':'; if ~isempty(varargin), I = varargin{1}; end
      P = obj.evalReferenceMap(points, 0, I); % nExnPxnW
      [nE, nP, nD] = size(P);
      P = reshape(P, nE*nP, nD); % (nE*nP)xnW
      switch nargin(F)
        case 1 % F(x)
          R = reshape(F(P), nE, nP, []);
        case 2 % F(x,U)
          if ~iscell(S.U); S.U = {S.U}; end
          for i = 1:numel(S.U)
            S.U{i} = reshape(S.U{i}, [], size(S.U{i},3));
          end
          R = reshape(F(P, S.U), nE, nP, []); % nExnPxnC
        case 3 % F(x,U,D)
          if ~iscell(S.U); S.U = {S.U}; end % nExnPxnC
          if ~iscell(S.dU); S.dU = {S.dU}; end % nExnPxnCxnD
          sz = size(S.dU{1});
          for i = 1:numel(S.U)
            S.U{i} = reshape(S.U{i}, [], sz(3));
            S.dU{i} = reshape(S.dU{i}, [], sz(3), sz(4));
          end
          R = reshape(F(P, S.U, S.dU), nE, nP, []); % nExnPxnCxnD
      end
    end
    function [R, RVec] = integrate(obj, func, quadRule, varargin)
      if ~isnumeric(func)
        func = obj.evalFunction(func, quadRule.points, [], varargin{:});
      end
      [~,~,trafo] = obj.evalTrafoInfo(quadRule.points,varargin{:});
      RVec = (func.*abs(trafo))*quadRule.weights;
      R =  sum(RVec);
    end
  end
  methods % refinement
    function uniformRefine(obj, N)
      for i = 1:N
        obj.nodes = obj.topology.uniformRefine()*obj.nodes;
      end
      obj.notifyObservers();
    end
    function uniformRefineFast(obj, N)
      for i = 1:N
        obj.nodes = obj.topology.uniformRefineFast()*obj.nodes;
      end
      obj.notifyObservers();
    end
  end
  methods % mesh operations
    function rotate(obj, alpha)
      if obj.dimW ~= 2, error('!Dimension must be 2!'); end
      A = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
      obj.applyLinearMap(A);
    end
    function scale(obj, a)
      A = diag(a)*eye(obj.dimW);
      obj.applyLinearMap(A);
    end
    function translate(obj, vec)
      obj.nodes = bsxfun(@plus, obj.nodes, vec(:)');
      obj.notifyObservers();
    end
    function applyLinearMap(obj, A)
      obj.nodes = (A*obj.nodes')';
      obj.notifyObservers();
    end
  end
  methods % mesh information
    function R = getMeasure(obj, dim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      ee = obj.topology.getEntity(dim); ee = ee(I,:);
      switch dim
        case 3
          v1 = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          v2 = obj.nodes(ee(:,3),:) - obj.nodes(ee(:,1),:);
          v3 = obj.nodes(ee(:,4+(size(ee,2)==8)),:) - obj.nodes(ee(:,1),:);
          R = ((v1(:,1).*v2(:,2).*v3(:,3) + v1(:,2).*v2(:,3).*v3(:,1)+v1(:,3).*v2(:,1).*v3(:,2)) ...
          - (v1(:,3).*v2(:,2).*v3(:,1)+v1(:,2).*v2(:,1).*v3(:,3)+v1(:,1).*v2(:,3).*v3(:,2)));
          if size(ee,2)==4
            R = R/6;  
          else
            if all(all(obj.nodes(ee(:,1),:) + v1+v2+v3 - obj.nodes(ee(:,8),:)>1e-12))
              warning('! Volume only valid for parallelepiped !');
            end
          end
        case 2
          v1 = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          v2 = obj.nodes(ee(:,3),:) - obj.nodes(ee(:,1),:);
          R = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1));
          if size(ee,2)==3
            R = R/2;
          else
            if all(all(obj.nodes(ee(:,1),:) + v1+v2 - obj.nodes(ee(:,4),:)>1e-12))
              warning('! Area only valid for parallelograms !');
            end
          end
        case 1
          v = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          R = sum(v.^2,2).^0.5;
       end
    end
    function R = getCenter(obj, dim, varargin) % [I]
      I = ':'; if nargin > 2, I = varargin{1}; end
      if dim == 0
        R = obj.nodes(I,:);
        return;
      end
      R = obj.topology.getEntity(dim); R = R(I,:);
      [nE, nV] = size(R);
      R = permute(mean(reshape(obj.nodes(R,:), nE, nV, []),2),[1 3 2]);
    end
    function R = getOuterNormal(obj, points)
      assert(obj.dimW==2);
      I = obj.isBoundary(); % bFace
      f2e = obj.topology.getFace2Elem(); f2e = f2e(I,:);
      R = obj.evalReferenceMap(points, 1, I); % nExnPx2
      R = reshape(([0 1;-1 0]*reshape(R,[],2)')',size(R));
      I = (f2e(:,1)==0);
      R(I,:,:) = -R(I,:,:); % nExnPx2
      R = bsxfun(@rdivide, R, sum(R.^2,3).^0.5); % nExnPx2
    end
    function R = getOuterNormal2(obj, element)
      fes = FESpace(obj, element);
      R = cell(1,obj.dimW);
      for k = 1:obj.dimW
        fc = FcDk(1,k,fes);
        fc.assemble();
        R{k} = fc.matrix();
      end
      R = cell2mat(R);
      R = bsxfun(@rdivide,R, sum(R.^2,2).^0.5);
    end
    function R = getDiam(obj)
      R = [min(obj.nodes); max(obj.nodes)];
    end
    function R = findEntity(obj, dim, varargin) % [loc]
      if nargin > 2
        loc = varargin{1};
        entity = obj.topology.getEntity(dim);
        R = any(reshape(loc(obj.nodes(entity,:)), size(entity,1), []), 2);
      else
        R = true(obj.topology.getNumber(dim), 1);
      end
    end
    function R = findEntityC(obj, dim, varargin) % [loc]
      if nargin > 2
        R = varargin{1}(obj.getCenter(dim)); % nE
      else
        R = true(obj.topology.getNumber(dim), 1);
      end
    end
    function R = isBoundary(obj, varargin) % [loc]
      R = obj.topology.isBoundary();
      if nargin > 1
        if ~isempty(varargin{1})
          I = varargin{1}(obj.getCenter(obj.topology.dimP-1, R));
          R = repmat(R, 1, size(I,2));
          R(R(:,1)>0,:) = I;
        else
          R = [];
          return
        end
      end
    end
    function R = getBoundary(obj, varargin) % [loc]
      R = obj.topology.getEntity('1', obj.isBoundary(varargin{:}));
    end
    function R = isSurface(obj, varargin) % [loc]
      if nargin < 2, R = obj.isBoundary(); return; end
      goodElem = obj.findEntity('0', varargin{:});
      E2F = obj.topology.getElem2Face();
      E2F = E2F(goodElem,:);
      uE2F = unique(E2F(:));
      I = hist(E2F(:),uE2F)==1;
      R = full(sparse(uE2F(I), 1, true, obj.topology.getNumber('1'), 1));
    end
    function R = getSurface(obj, varargin) % [loc]
      R = obj.topology.getEntity(1, obj.isSurface(varargin{:}));
    end
    function R = isBoundaryNode(obj, varargin) % [loc]
      R = unique(obj.getBoundary(varargin{:}));
      R = accumarray(R,1,[obj.topology.getNumber(0) 1])>0;
    end
    function R = getBoundaryNode(obj, varargin) % [loc]
      R = obj.topology.getEntity(0, obj.isBoundaryNode(varargin{:}));
    end
  end
  methods % mesh functions
    function R = getLocation(obj, dim, loc)
      R = loc(obj.getCenter(dim));
    end
    function showMeshFunction(obj, F, dim)
      if numel(F) ~= obj.topology.getNumber(dim)
        error('F is not mesh function of dimension dim');
      end
      if obj.topology.dimP~=2
        error('Mesh function supported for dim=2');
      end
      obj.show(); hold on;
      coord = obj.nodes;
      switch dim
        case 0
          plot3(coord(:,1),coord(:,2), F,'rx');
        case 1
          connect = obj.topology.connectivity{2,1}; nE = size(connect,1);
          X = [reshape(coord(connect,1),nE,[]) nan(nE,1)]';
          Y = [reshape(coord(connect,2),nE,[]) nan(nE,1)]';
          F = [repmat(F,1,size(connect,2)) nan(nE,1)]';
          plot3(X,Y,F,'r-')
        case 2
          connect = obj.topology.connectivity{3,1};
          nE = size(connect,1);
          X = reshape(coord(connect,1),nE,[]);
          X = [X X(:,1) nan(nE,1)]';
          Y = reshape(coord(connect,2),nE,[]);
          Y = [Y Y(:,1) nan(nE,1)]';
          F = [repmat(F,1,size(connect,2)+1) nan(nE,1)]';
          if ~obj.element.isSimplex()
            X([3 4],:) = X([4 3],:);
            Y([3 4],:) = Y([4 3],:);
            F([3 4],:) = F([4 3],:);
          end
          plot3(X,Y,F,'r-');
      end
      hold off
    end
  end
  methods % display
    function show(obj, varargin) % [type]
      v = VisualizerMesh.create(obj);
      v.show();
      if nargin < 2, return; end
      hold on;
      for i = 1:numel(varargin{1})
        v.showEntity(str2double(varargin{1}(i)));
      end
      hold off
      axis equal
    end
  end
  methods(Static = true)
    function R = getTopology(nodes, elem, dimP)
      switch size(elem, 2)
        case 2
          R = MeshTopologyInt(elem);
        case 3
          elem = MeshTopologyTri.renumber(nodes, elem);
          R = MeshTopologyTri(elem);
        case 4
          if dimP == 2
            R = MeshTopologyQuad(elem);
          else
            elem = MeshTopologyTet.renumber(nodes, elem);
            R = MeshTopologyTet(elem);
          end
        case 8
          R = MeshTopologyHex(elem);
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
    function [nodes, elem] = getMesh2D(grid, A, varargin)
      [nn, ee] = Mesh.getTensorProductMesh(grid, varargin{:});
      N = size(A);
      nodes = cell(prod(N),1);
      elem = cell(prod(N),1);
      cnt = -1;
      for i = 1:N(1)
        for j = 1:N(2)
          if ~A(i,j), continue; end
          cnt = cnt + 1;
          nodes{i+N(1)*(j-1)}(:,1) = nn(:,1) + (j-1)*grid{1}(end);
          nodes{i+N(1)*(j-1)}(:,2) = nn(:,2) + (N(1)-i)*grid{2}(end);
          elem{i+N(1)*(j-1)} = ee + cnt*size(nn,1);
        end
      end
      [nodes, ~, J] = uniquetol(cell2mat(nodes),1e-6,'byrows',1);
      elem = cell2mat(elem);
      elem = J(elem);
    end
    function R = transformTri2Quad(m)
      nSmooth = 10;
      %
      nN = m.topology.getNumber(0);
      nF = m.topology.getNumber(1);
      nE = m.topology.getNumber(2);
      [f2e, type] = m.topology.getFace2Elem();
      e2F = m.topology.getElem2Face();
      face = m.topology.getEntity(1);
      elem = m.topology.getEntity(2);
      % sort
      [~,I] = sort(m.getMeasure(1),'descend');
      % find pairs
      mem = false(nF,1); pairs = zeros(nF,3);
      cnt = 1;
      for k = 1:nF
        f = f2e(I(k),:);
        if prod(f)>0 && ~(mem(f(1)) || mem(f(2)))
          pairs(cnt,:) = [I(k) f];
          cnt = cnt + 1;
          mem(f,:) = true;
        end
      end
      pairs = pairs(1:cnt-1,:);
      single = find(~accumarray(reshape(pairs(:,[2 3]),[],1),1,[nE,1]));
      % refine
      Ep = cell(1,2);
      for cc = 1:2
        Ep{cc} = [elem(pairs(:,cc+1),:), nN+e2F(pairs(:,cc+1),:)];
        I = type(pairs(:,1),cc)==1; Ep{cc}(I,:) = Ep{cc}(I,[3 1 2 6 4 5]);
        I = type(pairs(:,1),cc)==2; Ep{cc}(I,:) = Ep{cc}(I,[1 2 3 4 5 6]);
        I = type(pairs(:,1),cc)==3; Ep{cc}(I,:) = Ep{cc}(I,[2 3 1 5 6 4]);
      end
      Ep = cell2mat(Ep);
      Es = [elem(single,:), nN+e2F(single,:), nN+nF+(1:numel(single))'];
      E = [Ep(:,[1 4 6 5]);Ep(:,[4 2 5 12]);Ep(:,[6 5 3 10]);Ep(:,[11 12 10 7]); ...
           Es(:,[1 4 6 7]);Es(:,[4 2 7 5]);Es(:,[7 5 6 3])];
      % new nodes
      P = [speye(nN); sparse(repmat((1:nF)',1,2), face, 0.5, nF, nN); ...
                      sparse(repmat((1:numel(single))',1,3), elem(single,:), 1/3, numel(single), nN)];
      NN = P*m.nodes;
      %
      R = Mesh(NN, E);
      for k = 1:nSmooth
        zNN = [0 0; R.nodes];
        n2N = R.topology.getNodePatch(0);
        xij = reshape(zNN(n2N+1,:),size(NN,1),[],2);
        %dij = sum((xij - permute(zNN(2:end,:),[1 3 2])).^2,3).^0.5; dij(n2N==0) = 0;
        NN = permute(sum(xij,2),[1 3 2])./sum(abs(xij(:,:,1))>0,2);
        R.nodes(~R.isBoundaryNode(),:) = NN(~R.isBoundaryNode(),:);
      end
    end
  end
end