classdef FESpace < SOFEClass
  properties
    mesh, element, quadRule
    fixB, shift, freeDoFs
    cache, isCaching = true;
    observers
  end
  methods % constructor, observer and caching.
    function obj = FESpace(mesh, element, varargin) % [fixB, shift]
      if mesh.element.isSimplex ~= element.isSimplex || ...
         mesh.element.dimension ~= element.dimension
        warning('! Mesh and Element are not compatible, continue? !');
        keyboard
      end
      obj.mesh = mesh;
      obj.element = element;
      obj.quadRule = obj.mesh.topology.getQuadRule(max(2*(obj.element.order),1));
      obj.mesh.topology.register(obj);
      obj.observers = {};
      obj.resetCache();
      obj.fixB = @(x)false(size(x,1),1);
      obj.shift = [];
      if nargin > 2
        obj.fixB = varargin{1};
        if nargin > 3
          obj.shift = varargin{2};
        end
      end
      obj.freeDoFs = ~obj.getBoundaryDoFs(obj.fixB);
    end
    function resetCache(obj)
      nBlock = obj.mesh.nBlock;
      obj.cache.refMaps = cell(nBlock, 1);
      obj.cache.DPhi = cell(nBlock, 1);
      obj.cache.DPhiInv = cell(nBlock, 1);
      obj.cache.jac = cell(nBlock, 1);
      obj.cache.basis = cell(nBlock, 1);
      for k = 1:nBlock
        obj.cache.refMaps{k} = cell(obj.element.dimension+1,1);
        obj.cache.DPhi{k} = cell(obj.element.dimension+1,1);
        obj.cache.DPhiInv{k} = cell(obj.element.dimension+1,1);
        obj.cache.jac{k} = cell(obj.element.dimension+1,1);
        obj.cache.basis{k} = cell(obj.element.dimension+1,3);
      end
      obj.cache.dM = [];
    end
    function R = getBlock(obj, codim, varargin)
      R = obj.mesh.getBlock(codim, varargin{:});
    end
    function notify(obj)
      obj.resetCache();
      obj.freeDoFs = ~obj.getBoundaryDoFs(obj.fixB);
      obj.notifyObservers();
    end
    function register(obj, observer)
      obj.observers = [obj.observers, {observer}];
    end
    function notifyObservers(obj)
      for i = 1:numel(obj.observers)
        obj.observers{i}.notify();
      end
    end
  end
  methods % quadrature.
    function setQuadRule(obj, quadRule)
      obj.quadRule = quadRule;
      obj.resetCache();
    end
    function [Rp, Rw] = getQuadData(obj, codim)
      Rp = obj.quadRule{codim+1}.points;
      Rw = obj.quadRule{codim+1}.weights;
    end
  end
  methods % evaluation.
    function R = evalFunction(obj, F, points, codim, U, varargin) % [k or I]
      if isempty(points)
        k = varargin{1};
        idx = obj.getBlock(codim, k);
        varargin{1} = (idx(1):idx(2))';
        points = obj.getQuadData(codim);
      end
      R = obj.mesh.evalFunction(F, points, U, varargin{:});
    end
    function R = evalReferenceMap(obj, points, codim, varargin) % [k or I]
      if isempty(points)
        if nargin < 4
          R = [];
          for k = 1:obj.mesh.nBlock
            R = [R; obj.evalReferenceMap(points, codim, k)];
          end
          return
        end
        k = varargin{1};
        if ~isempty(obj.cache.refMaps{k}{codim+1}) && obj.isCaching
          R = obj.cache.refMaps{k}{codim+1};
          return
        else
          idx = obj.getBlock(codim, k);
          R = obj.mesh.evalReferenceMap(obj.getQuadData(codim), 0, (idx(1):idx(2))');
          if isempty(points) && obj.isCaching
            obj.cache.refMaps{k}{codim+1} = R;
          end
        end
      else
        R = obj.mesh.evalReferenceMap(points, 0, varargin{:});
      end
    end
    function [R, invR, jacR] = evalTrafoInfo(obj, points, codim, varargin) % [k or I]
      if isempty(points)
        k = varargin{1};
        if ~isempty(obj.cache.DPhi{k}{codim+1}) && obj.isCaching
          R = obj.cache.DPhi{k}{codim+1};
          invR = obj.cache.DPhiInv{k}{codim+1};
          jacR = obj.cache.jac{k}{codim+1};
          return
        else
          idx = obj.getBlock(codim, k);
          [R, invR, jacR] = obj.mesh.evalTrafoInfo(obj.getQuadData(codim), (idx(1):idx(2))');
          if obj.isCaching
            obj.cache.DPhi{k}{codim+1} = R;
            obj.cache.DPhiInv{k}{codim+1} = invR;
            obj.cache.jac{k}{codim+1} = jacR;
          end
        end
      else
        [R, invR, jacR] = obj.mesh.evalTrafoInfo(points, varargin{:});
      end
    end
    function R = evalDoFVector(obj, U, points, codim, order, varargin) % [k or I]
      assert(numel(U)==obj.getNDoF(), 'First argument must be DoFVector!');
      if iscell(points) % global evaluation
        if numel(points) == 1
          points = obj.mesh.evalInversReferenceMap(points{1});
        end
        isValid = points{2}>0;
        points{1} = points{1}(isValid,:); points{2} = points{2}(isValid);
        basis = obj.evalGlobalBasis(points, 0, order, points{2}); % [1/nE]xnB[xnP]xnCx[nD]       
        dMap = abs(obj.getDoFMap(0, points{2})).'; % nExnB
        value = zeros(size(dMap)); % nExnB
        I = dMap > 0; value(I) = U(dMap(I));
        value = sum(bsxfun(@times, permute(value,[1 3:4 2]), ...
                                   permute(basis,[1 3:4 2])), 4); % nExnCx[nD]
        R = nan(numel(isValid), size(basis, 3), size(basis, 4)); % nExnC[xnD]
        R(isValid,:,:) = value; % nExnC[xnD]
      else % local evaluation
        if ~isempty(points)
          codim = obj.element.dimension - size(points,2);
          if nargin >5, idx = varargin{:}; else idx = ':'; end
        else
          idx = obj.getBlock(codim, varargin{1});
          idx = (idx(1):idx(2))';
        end
        assert(codim==0 || order==0, '! Derivatives for traces not supported !');
        basis = obj.evalGlobalBasis(points, codim, order, varargin{:}); % [1/nE]xnB[xnP]xnCx[nD]
        dMap = abs(obj.getDoFMap(codim, idx)).'; % nExnB
        R = zeros(size(dMap)); % nExnB
        I = dMap > 0; R(I) = U(dMap(I));
        R = sum(bsxfun(@times, permute(R,[1 3:5 2]), ...
                               permute(basis,[1 3:5 2])), 5); % nExnPxnCx[nD]
      end
    end
    function R = evalGlobalBasis(obj, points, codim, order, varargin) % [k or I]
      if isempty(points)
        if nargin > 4 && isnumeric(varargin{1}) % not ':'
          k = varargin{1};
          if ~isempty(obj.cache.basis{k}{codim+1, order+1}) && obj.isCaching
            R = obj.cache.basis{k}{codim+1, order+1};
          else
            R = obj.computeGlobalBasis(points, codim, order, k);
            if obj.isCaching
              obj.cache.basis{k}{codim+1, order+1} = R;
            end
          end
        else % block evaluation
          nBlock = size(obj.getBlock(codim),2);
          R = cell(nBlock,1);
          for k = 1:nBlock;
            if ~isempty(obj.cache.basis{k}{codim+1, order+1}) && obj.isCaching
              R{k} = obj.cache.basis{k}{codim+1, order+1};
            else
              R{k} = obj.computeGlobalBasis(points, codim, order, k);
              if obj.isCaching
                obj.cache.basis{k}{codim+1, order+1} = R{k};
              end
            end
          end
          R = cell2mat(R);
        end
      else
        if nargin > 4 || iscell(points)
          R = obj.computeGlobalBasis(points, codim, order, varargin{:});
        else % block evaluation
          idx = obj.getBlock(obj.element.dimension-size(points,2));
          R = cell(size(idx,2),1);
          for k = 1:size(idx,2)
            R{k} = obj.computeGlobalBasis(points, codim, order, (idx(1,k):idx(2,k))');
          end
          R = cell2mat(R);
        end
      end
    end
    function R = computeGlobalBasis(obj, points, codim, order, varargin) % [k or I]
      [trafo{1},trafo{2},trafo{3}] = obj.evalTrafoInfo(points, codim, varargin{:}); % nExnP[xnWxnD]
      if iscell(points)
        basis = obj.element.evalBasis(points{1}, order); % nBxnPxnC[xnD]
        trafo{1} = permute(trafo{1}, [4 5 1 2 3]); % 1x1xnPxnWxnD
        trafo{2} = permute(trafo{2}, [4 5 1 2 3]); % 1x1xnPxnDxnW
        trafo{3} = permute(trafo{3}, [2 3 1]); % 1x1xnP
        I = points{2}; codim = obj.element.dimension - size(points{1},2);
      else
        if isempty(points)
          idx = obj.getBlock(codim, varargin{1}); I = (idx(1):idx(2))';
          basis = obj.element.evalBasis(obj.getQuadData(codim), order); % nBxnPxnC[xnD]
        else
          if nargin > 4, I = varargin{1}; else I = ':'; end
          basis = obj.element.evalBasis(points, order); % nBxnPxnC[xnD]
          codim = obj.element.dimension - size(points,2);
        end
        trafo{1} = permute(trafo{1}, [1 5 2 3 4]); % nEx1xnPxnWxnD
        trafo{2} = permute(trafo{2}, [1 5 2 3 4]); % nEx1xnPxnDxnW
        trafo{3} = permute(trafo{3}, [1 3 2]); % nEx1xnP
      end
      basis = permute(basis, [5 1 2 3 4]); % 1xnBxnPxnC[xnD]
      switch obj.element.conformity
        case 'L2'
          R = basis; % 1xnBxnPxnC[xnD]
        case 'H1'
          switch order
            case 0
              R = basis; % 1xnBxnPxnC[xnD]
            case 1
              R = sum(bsxfun(@times, permute(basis,    [1 2 3 4 6 5]), ...
                                     permute(trafo{2}, [1 2 3 6 5 4])), 6); % nExnBxnPxnCxnW
          end
        case 'HDiv'
          switch order
            case 0
              R = sum(bsxfun(@times, trafo{1}, permute(basis, [1 2 3 5 4])),5); % nExnBxnPxnW
              R = bsxfun(@ldivide, trafo{3}, R); % nExnBxnPxnW
            case 1
              R = sum(bsxfun(@times, permute(basis, [1 2 3 4 6 5]), ...
                                     permute(trafo{2}, [1 2 3 6 5 4])), 6); % nExnBxnPxnCxnW
              R = sum(bsxfun(@times, permute(trafo{1}, [1 2 3 4 6 5]), ...
                                     permute(R, [1 2 3 6 5 4])), 6); % nExnBxnPxnWxnW
              R = bsxfun(@ldivide, trafo{3}, R); % nExnBxnPxnWxnW
          end
        case 'HRot'
          switch order
            case 0
              R = sum(bsxfun(@times, permute(trafo{2}, [1 2 3 5 4]), ...
                                     permute(basis, [1 2 3 5 4])), 5); % nExnBxnPxnW
            case 1
              R = sum(bsxfun(@times, permute(basis, [1 2 3 4 6 5]), ...
                                     permute(trafo{2}, [1 2 3 6 5 4])), 6); % nExnBxnPxnCxnW
              R = sum(bsxfun(@times, permute(trafo{2}, [1 2 3 5 6 4]), ...
                                     permute(R, [1 2 3 6 5 4])), 6); % nExnBxnPxnWxnW
          end
      end
      if iscell(points)
        R = permute(R,[3 2 4 5 1]); % nExnBx[...]
      end
      dMap = obj.getDoFMap(codim, I).'; % nExnB
      R = bsxfun(@times, sign(dMap), R); % nExnBxnPxnC[xnD]
    end
  end
  methods % DoFManager.
    function [R, nDoF] = getDoFMap(obj, codim, varargin) % [I]
      if ~isempty(obj.cache.dM)
        R = obj.cache.dM.doFArrays;
        nDoF = obj.cache.dM.nDoF;
      else
        [R, nDoF] = obj.computeDoFMaps();
        obj.cache.dM.doFArrays = R;
        obj.cache.dM.nDoF = nDoF;
      end
      R = R{codim+1}; % nBxnE
      if nargin > 2
        R = R(:,varargin{1});
      end
    end
    function [R, nDoF] = computeDoFMaps(obj)
      fprintf('Compute DoF maps ... ');
      doFs = cell(obj.element.dimension+1,1);
      R = cell(size(doFs));
      nDoF = 0;
      for dim = 0:numel(R)-1 % assemble doFMaps
        R{dim+1} = cell(dim+1,1);
        N = obj.mesh.topology.getNumber(dim) * prod(obj.element.doFTuple(:,dim+1));
        doFs{dim+1} = reshape(nDoF+(1:N), [], obj.mesh.topology.getNumber(dim)); % nDoFLocxnE
        nDoF = nDoF + N;
        for d = 0:dim
          C = obj.mesh.topology.connectivity{dim+1, d+1}; % nExnESub
          R{dim+1}{d+1} = reshape(obj.orient(doFs{d+1}(:,C'), dim, d), [], size(C,1)); % {nD}{nD}(nDoFLoc*nESub)xnE
        end
        R{dim+1} = cell2mat(R{dim+1}); % {nD}(nBxnE)
      end
      R = R(end:-1:1); % dim+1 --> codim+1
      fprintf('DONE\n');
    end
    function R = orient(obj, R, dim, d)
      if dim < 2, return; end
      orient = permute(obj.mesh.topology.getOrientation(dim, d),[2 1 3]); % nESubxnExnO
      orient = reshape(orient, [], size(orient, 3)); % (nESub*nE)xnO
      tmp = zeros(size(R)); % nDoFLocx(nESub*nE)
      if d == 1 && dim > d % EDGE
        for k = -1:2:1
          I = obj.element.getDoFEnum(1,k); % nDoFLoc
          if ~isempty(I)
            iO = orient == k; % nESub*nE
            tmp(abs(I),iO) = bsxfun(@times, sign(I(:)), R(:, iO)); % nDoFLocx(nESub*nE)
          end
        end
        R = tmp; % nDoFLocx(nESub*nE)
      end
      if d == 2 && dim > d % FACE
        if obj.element.isSimplex()
          for s = 3:-1:1
            for k = -1:2:1
              I = obj.element.getDoFEnum(2,s,k); % nDoFLoc
              if ~isempty(I)
                iO = all(bsxfun(@eq, orient, [s k]),2); % nESub*nE
                tmp(abs(I),iO) = bsxfun(@times, sign(I), R(:, iO)); % nDoFLocx(nESub*nE)
              end
            end
          end
        else
          for s1 = -1:2:1
            for s2 = -1:2:1
              for k = -1:2:1
                I = obj.element.getDoFEnum(2,s1,s2,k); % nDoFLoc
                if ~isempty(I)
                  iO = all(bsxfun(@eq, orient, [s1 s2 k]),2); % (nESub*nE)
                  tmp(abs(I),iO) = bsxfun(@times, sign(I(:)), R(:, iO)); % nDoFLocx(nESub*nE)
                end
              end
            end
          end
        end
        R = tmp; % nDoFLocx(nESub*nE)
      end
      if dim == obj.element.dimension && d == dim-1 && strcmp(obj.element.conformity, 'HDiv')
        orient = reshape(obj.mesh.topology.getNormalOrientation(dim, d).',1,[]); % 1x(nESub*nE)
        R = bsxfun(@times, R, orient); % nDoFLocx(nESub*nE)
      end
    end
    function nDoF = getNDoF(obj)
      [~, nDoF] = obj.getDoFMap(0);
    end
    function R = extractDoFs(obj, codim, I)
      [dMap, nDoF] = obj.getDoFMap(codim);
      %
      nC = size(I,2);
      doFs = cell(nC,1);
      for d = 1:nC
        doFs{d} = reshape(abs(dMap(d:nC:end, I(:,d))), [], 1);
      end
      doFs = cell2mat(doFs);
      %
      R = accumarray(setdiff(unique(doFs(:)), 0), 1, [nDoF 1])>0;
    end
    function R = getBoundaryDoFs(obj, varargin) % [loc]
      R = obj.extractDoFs(1, obj.mesh.topology.isBoundary(varargin{:}));
    end
  end
  methods % interpolation.
    function R = getShiftVector(obj, varargin) % [time]
      R = zeros(obj.getNDoF,1);
      if ~isempty(obj.shift)
        if nargin(obj.shift) > 1
          if nargin > 1
            func = @(x)obj.shift(x, varargin{1});
          else
            func = @(x)obj.shift(x, 0.0);
          end
        else
          func = obj.shift;
        end
        if obj.mesh.element.dimension == 1
          R = obj.getWeakInterpolation(func, 0); % nDoFx1
        else
          R = obj.getWeakInterpolation(func, 1, obj.mesh.topology.isBoundary()); % nDoFx1
        end
      end
    end
    function R = getWeakInterpolation(obj, f, codim, varargin)
      [P, weights] = obj.getQuadData(codim); % nPx1
      basis = obj.evalGlobalBasis(P, [], 0, varargin{:}); % nExnBxnPxnC
      F = permute(obj.evalFunction(f, P, [], [], varargin{:}), [1 4 2 3]); % nEx1xnPxnC
      [~,~,jac] = obj.mesh.evalTrafoInfo(P, varargin{:}); % nExnP
      dX = bsxfun(@times, abs(jac), weights'); % nExnP
      lhs = sum(bsxfun(@times, permute(basis,[1 2 5 3 4]), ...
                               permute(basis,[1 5 2 3 4])), 5); % nExnBxnBxnP
      lhs = permute(sum(bsxfun(@times, lhs, permute(dX, [1 3 4 2])), 4),[2 3 1]); % nBxnBxnE
      rhs = sum(bsxfun(@times, basis, F), 4); % nExnBxnP
      rhs = sum(bsxfun(@times, rhs, permute(dX, [1 3 2])),3).'; % nBxnE
      %
      dMap = reshape(abs(obj.getDoFMap(codim, varargin{:})),[],1); % nB*nE
      iZ = (dMap~=0); dMap = abs(dMap(:));
      if any(~iZ(:))
        dMap = dMap(iZ);
        rhs = rhs(iZ);
        [nB,~, nE] = size(lhs);
        lhs = reshape(lhs,nB,[]);
        lhs = reshape(lhs(:,iZ),nB,[],nE);
        [~,nBNew,~] = size(lhs);
        lhs = reshape(permute(lhs,[2 1 3]),nBNew,(nB*nE)); % nBNew*(nBxnE)
        lhs = reshape(lhs(:,iZ),nBNew,nBNew,nE);
      end
      lhs = blockify(lhs); % (nB*nE)x(nB*nE)
      rhs = lhs\rhs(:); % nB*nE
      R = zeros(obj.getNDoF(),1); % nDoFx1
      rhs = accumarray(dMap,rhs)./accumarray(dMap,1);
      I = unique(dMap);
      R(I) = rhs(I);
    end
    function R = getStrongInterpolation(obj, f, codim, varargin)
      P = obj.getLocalEquiPoints(obj.element.dimension - codim);
      basis = obj.evalGlobalBasis(P, [], 0, varargin{:}); % nExnBxnPxnC
      lhs = permute(reshape(basis,size(basis,1),size(basis,2),[]), [3 2 1]); % (nP*nC)xnBxnE
      rhs = permute(obj.evalFunction(f, P, [], [], varargin{:}), [2 3 1]); % nP*nC*nE
      %
      dMap = reshape(abs(obj.getDoFMap(codim, varargin{:})),[],1); % nB*nE
      iZ = (dMap~=0); dMap = abs(dMap(:));
      if any(~iZ(:))
        dMap = dMap(iZ);
        rhs = rhs(iZ);
        [nB,~, nE] = size(lhs);
        lhs = reshape(lhs,nB,[]);
        lhs = reshape(lhs(:,iZ),nB,[],nE);
        [~,nBNew,~] = size(lhs);
        lhs = reshape(permute(lhs,[2 1 3]),nBNew,(nB*nE)); % nBNew*(nBxnE)
        lhs = reshape(lhs(:,iZ),nBNew,nBNew,nE);
      end
      lhs = blockify(lhs); % (nB*nE)x(nB*nE)
      rhs = lhs\rhs(:); % nB*nE
      R = zeros(obj.getNDoF(),1); % nDoFx1
      rhs = accumarray(dMap,rhs)./accumarray(dMap,1);
      I = unique(dMap);
      R(I) = rhs(I);
    end
    function R = getGlobalInterpolation(obj, f)
      mass = Op_data_Id_Id(1, 0, obj); mass.assemble();
      l2 = Fc_Data_Id(f, obj,0); l2.assemble();
      R = mass.matrix \ l2.vector;
    end
    function R = getLocalEquiPoints(obj, varargin) % [nD]
      if nargin > 1, nD = varargin{1}; else nD = obj.element.dimension; end
      points = linspace(0,1,obj.element.order+1)';
      R = points;
      for i = 2:nD
        R = [kron(ones(length(points),1),R), kron(points,ones(length(points)^(i-1),1))];
      end
      if obj.element.isSimplex()
        R = R((sum(R,2)<=1),:);
      end
    end
  end
end
