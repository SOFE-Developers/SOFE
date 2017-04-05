classdef FESpace < SOFEClass
  properties
    mesh, element, quadRule
    fixB, shift, freeDoFs
    cache, isCaching = true;
    observers
    nBlock = 1;
  end
  methods % constructor, observer and caching.
    function obj = FESpace(mesh, element, varargin) % [fixB, shift]
      try
        if mesh.element.isSimplex ~= element.isSimplex || ...
           mesh.element.dimension ~= element.dimension
          warning('! Mesh and Element are not compatible, continue? !');
          keyboard
        end
      end
      obj.mesh = mesh;
      obj.element = element;
      obj.quadRule = obj.mesh.topology.getQuadRule(max(2*(obj.element.order),1));
      obj.setBlocking();
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
    end
    function resetCache(obj, varargin)
      nBlock = max(obj.nBlock);
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
      if nargin < 2
        obj.cache.dM = [];
        obj.freeDoFs = [];
      end
    end
    function setBlocking(obj)
      nC = obj.element.getNC(); nD = obj.element.dimension;
      nB = obj.element.nB(end); nQ = numel(obj.quadRule{1}.weights);
      elPerBlock = SOFEClass.getElementsPerBlock(nB, nQ, nC, nD);
      obj.nBlock = ceil(obj.mesh.topology.getNumber(nD)/elPerBlock);
      nB = obj.element.nB(end-1); nQ = numel(obj.quadRule{2}.weights);
      elPerBlock = SOFEClass.getElementsPerBlock(nB, nQ, nC, nD);
      obj.nBlock(2) = ceil(obj.mesh.topology.getNumber(nD-1)/elPerBlock);
    end
    function R = getBlock_(obj, codim, varargin)
      R = obj.mesh.getBlock(codim, varargin{:});
    end
    function R = getBlock(obj, codim, varargin) % [k]
      nE = obj.mesh.topology.getNumber(obj.mesh.topology.dimP - codim);
      if obj.nBlock(codim+1)>nE, error('!Number of blocks exceeds number of elements!'); end
      R = unique(floor(linspace(0,nE,obj.nBlock(codim+1)+1)));
      R = [R(1:end-1)+1; R(2:end)];
      if nargin > 2
        R = (R(1,varargin{:}):R(2,varargin{:}))';
      else
        R = (1:R(end))';
      end
    end
    function notify(obj)
      obj.setBlocking();
      obj.resetCache();
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
      obj.resetCache('Do not reset DoFMaps');
    end
    function [Rp, Rw] = getQuadData(obj, codim)
      Rp = obj.quadRule{codim+1}.points;
      Rw = obj.quadRule{codim+1}.weights;
    end
  end
  methods % evaluation.
    function R = evalFunction(obj, F, points, codim, U, varargin) % [{k} or I]
      block = false;
      if nargin > 5 && iscell(varargin{1})
        block = true; k = varargin{1};
      end
      if ~isempty(points)
        codim = obj.element.dimension-size(points,2);
      end
      if block
        varargin{1} = obj.getBlock(codim, k{1});
      end
      if isempty(points)
        points = obj.getQuadData(codim);
      end
      R = obj.mesh.evalFunction(F, points, U, varargin{:});
    end
    function R = evalReferenceMap(obj, points, codim, varargin) % [{k} or I]
      block = false;
      if nargin > 3 && iscell(varargin{1})
        block = true; k = varargin{1};
      end
      cache =  obj.isCaching && block && isempty(points);
      %
      if cache && ~isempty(obj.cache.refMaps{k{1}}{codim+1})
        R = obj.cache.refMaps{k{1}}{codim+1};
      else
        if ~isempty(points)
          codim = obj.element.dimension-size(points,2);
        end
        if nargin < 4 || (nargin == 4 && ischar(varargin{1}) && strcmp(varargin,':'))
          nBlock = obj.nBlock(codim+1); R = cell(nBlock,1);
          for k = 1:nBlock
            R{k} = obj.evalReferenceMap(points, codim, {k});
          end
          try
            R = cell2mat(R);
          catch
            R = padcell2mat(R);
          end
          [~,I] = sort(obj.getBlock(codim));
          R = R(I,:,:,:);
        else
          if block
            varargin{1} = obj.getBlock(codim, k{1});
          end
          if isempty(points)
            points = obj.getQuadData(codim);
          end
          R = obj.mesh.evalReferenceMap(points, 0, varargin{:});
        end
        if cache
          obj.cache.refMaps{k{1}}{codim+1} = R;
        end
      end
    end
    function [R, invR, jacR] = evalTrafoInfo(obj, points, codim, varargin) % [{k} or I]
      block = false;
      if nargin > 3 && iscell(varargin{1})
        block = true; k = varargin{1};
      end
      cache =  obj.isCaching && block && isempty(points);
      %
      if cache && ~isempty(obj.cache.DPhi{k{1}}{codim+1})
        R = obj.cache.DPhi{k{1}}{codim+1};
        invR = obj.cache.DPhiInv{k{1}}{codim+1};
        jacR = obj.cache.jac{k{1}}{codim+1};
      else
        if ~isempty(points)
          codim = obj.element.dimension-size(points,2);
        end
        if block
          varargin{1} = obj.getBlock(codim, k{1});
        end
        if isempty(points)
          points = obj.getQuadData(codim);
        end
        [R, invR, jacR] = obj.mesh.evalTrafoInfo(points, varargin{:});
        if cache
          obj.cache.DPhiInv{k{1}}{codim+1} = invR;
          obj.cache.jac{k{1}}{codim+1} = jacR;
        end
      end
    end
    function R = computeGlobalBasis(obj, points, codim, order, varargin) % [{k} or I]
      block = false;
      if nargin > 4 && iscell(varargin{1})
        block = true; k = varargin{1};
      end
      if iscell(points)
        varargin{1} = points{2}; codim = obj.element.dimension - size(points{1},2);
        basis = obj.element.evalBasis(points{1}, order); % nBxnPxnC[xnD]
        pVec{2} = [2 3 1]; pVec{1} = [4 5 1 2 3];
      else
        if ~isempty(points)
          codim = obj.element.dimension-size(points,2);
        end
        if block
          varargin{1} = obj.getBlock(codim, k{1});
        end
        if isempty(points)
          points = obj.getQuadData(codim);
        end
        basis = obj.element.evalBasis(points, order); % nBxnPxnC[xnD]
        pVec{2} = [1 3 2]; pVec{1} = [1 5 2 3 4];
      end
      if ~(order==0 && (strcmp(obj.element.conformity, 'H1') || strcmp(obj.element.conformity, 'L2')))
        [trafo{1},trafo{2},trafo{3}] = obj.evalTrafoInfo(points, codim, varargin{:}); % nExnP[xnWxnD]
        trafo{1} = permute(trafo{1}, pVec{1}); % nEx1xnPxnWxnD
        trafo{2} = permute(trafo{2}, pVec{1}); % nEx1xnPxnDxnW
        trafo{3} = permute(trafo{3}, pVec{2}); % nEx1xnP
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
      dMap = obj.getDoFMap(codim, varargin{:}).'; % nExnB
      R = bsxfun(@times, sign(dMap), R); % nExnBxnPxnC[xnD]
    end
    function R = evalGlobalBasis(obj, points, codim, order, varargin) % [{k} or I]
      block = false;
      if nargin > 4 && iscell(varargin{1})
        block = true; k = varargin{1};
      end
      cache =  obj.isCaching && block && isempty(points);
      if cache && ~isempty(obj.cache.basis{k{1}}{codim+1, order+1})
        R = obj.cache.basis{k{1}}{codim+1, order+1};
      else
        if ~isempty(points)
          codim = obj.element.dimension-size(points,2);
        end
        if nargin < 5 || (nargin == 5 && ischar(varargin{1}) && strcmp(varargin,':'))
          nBlock = obj.nBlock(codim+1); R = cell(nBlock,1);
          for k = 1:nBlock
            R{k} = obj.evalGlobalBasis(points, codim, order, {k});
          end
          try
            R = cell2mat(R);
          catch
            R = padcell2mat(R);
          end
          [~,I] = sort(obj.getBlock(codim));
          R = R(I,:,:,:);
        else
          if block
            varargin{1} = obj.getBlock(codim, k{1});
          end
          if isempty(points)
            points = obj.getQuadData(codim);
          end
          R = obj.computeGlobalBasis(points, [], order, varargin{:});
        end
        if cache
          obj.cache.basis{k{1}}{codim+1, order+1} = R;
        end
      end
    end
    function R = evalDoFVector(obj, U, points, codim, order, varargin) % [{k} or I]
      assert(numel(U)==obj.getNDoF(), 'First argument must be DoFVector!');
      if iscell(points)
        R = obj.evalDoFVectorGlobal(U, points, codim, order);
      else
        R = obj.evalDoFVectorLocal(U, points, codim, order, varargin{:});
      end
    end
    function R = evalDoFVectorGlobal(obj, U, points, codim, order)
      if numel(points) == 1
        points = obj.mesh.evalInversReferenceMap(points{1});
      end
      isValid = points{2}>0;
      points{1} = points{1}(isValid,:); points{2} = points{2}(isValid);
      basis = obj.evalGlobalBasis(points, 0, order, []); % [1/nE]xnB[xnP]xnCx[nD]       
      dMap = abs(obj.getDoFMap(0, points{2})).'; % nExnB
      value = zeros(size(dMap)); % nExnB
      iZ = dMap > 0; value(iZ) = U(dMap(iZ));
      value = sum(bsxfun(@times, permute(value,[1 3:4 2]), ...
                                 permute(basis,[1 3:4 2])), 4); % nExnCx[nD]
      R = nan(numel(isValid), size(basis, 3), size(basis, 4)); % nExnC[xnD]
      R(isValid,:,:) = value; % nExnC[xnD]
    end
    function R = evalDoFVectorLocal(obj, U, points, codim, order, varargin) % [{k} or I]
        if ~isempty(points)
          codim = obj.element.dimension-size(points,2);
        end
        if nargin < 6 || (nargin == 6 && ischar(varargin{1}) && strcmp(varargin,':'))
        nBlock = obj.nBlock(codim+1); R = cell(nBlock,1); s = 0;
        for k = 1:nBlock
          R{k} = obj.evalDoFVector(U, points, codim, order, {k}); % nExnPxnCxnD
          if nBlock>1 
            fprintf(repmat('\b',1,length(s)));
            s = sprintf('progress evalDoFVector: %d / %d', k, nBlock); fprintf(s);
          end
        end
        if nBlock>1, fprintf('\n'); end
        try
          R = cell2mat(R);
        catch
          R = padcell2mat(R);
        end
        if ~isempty(points)
          codim = obj.element.dimension-size(points,2);
        end
        [~,I] = sort(obj.getBlock(codim));
        R = R(I,:,:,:);
      else
        block = false;
        if nargin > 5 && iscell(varargin{1})
          block = true; k = varargin{1};
        end
        if block
          varargin{1} = obj.getBlock(codim, k{1});
        end
        if isempty(points)
          points = obj.getQuadData(codim);
        end
        assert(codim==0 || order==0, '! Derivatives for traces not supported !');
        basis = obj.evalGlobalBasis(points, [], order, varargin{:}); % [1/nE]xnB[xnP]xnCx[nD]
        dMap = abs(obj.getDoFMap(codim, varargin{:})).'; % nExnB
        R = zeros(size(dMap)); I = dMap > 0; R(I) = U(dMap(I)); % nExnB
        R = sum(bsxfun(@times,permute(R,[1 3:5 2]),permute(basis,[1 3:5 2])),5); % nExnPxnCx[nD]
      end
    end
  end
  methods % DoFManager.
    function [R, nDoF] = getDoFMap(obj, codim, varargin) % [{k} or I]
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
        I = varargin{1};
        if iscell(I)
          I = obj.getBlock(codim, I{1});
        end
        R = R(:,I);
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
    function R = orient(obj, R, dim, d, varargin) % [I]
      if dim < 2, return; end
      orient = permute(obj.mesh.topology.getOrientation(dim, d, varargin{:}),[2 1 3]); % nESubxnExnO
      orient = reshape(orient, [], size(orient, 3)); % (nESub*nE)xnO
      tmp = zeros(size(R)); % nDoFLocx(nESub*nE)
      if d == 1 && dim > d % EDGE
        for k = -1:2:1
          I = obj.element.getDoFEnum(1,k); % nDoFLoc
          if ~isempty(I)
            iO = orient == k; % nESub*nE
            tmp(abs(I),iO) = bsxfun(@times, sign(I(:)), R(1:numel(I), iO)); % nDoFLocx(nESub*nE)
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
                iO = orient==(k*s); % nESub*nE
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
      tmp = setdiff(unique(doFs(:)), 0);
      R = accumarray(tmp(:), 1, [nDoF 1])>0;
    end
    function R = getBoundaryDoFs(obj, varargin) % [loc]
      R = obj.extractDoFs(1, obj.mesh.topology.isBoundary(varargin{:}));
    end
    function R = getFreeDoFs(obj)
      if ~isempty(obj.freeDoFs)
        R = obj.freeDoFs;
      else
        R = ~obj.getBoundaryDoFs(obj.fixB);
        obj.freeDoFs = R;
      end
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
          R = obj.getInterpolation(func, 0); % nDoFx1
        else
          R = obj.getInterpolation(func, 1, obj.mesh.topology.isBoundary()); % nDoFx1
        end
      end
    end
    function R = getL2Projection(obj, f)
      mass = Op_data_Id_Id(1, 0, obj); mass.assemble();
      l2 = Fc_Data_Id(f, obj,0); l2.assemble();
      R = mass.matrix \ l2.vector;
    end
    function R = getInterpolation(obj, f, codim, varargin) % [{k} or I]
      if nargin < 4 || (nargin == 6 && ischar(varargin{1}) && strcmp(varargin,':'))
        nBlock = obj.nBlock(codim+1); R = cell(1,nBlock); s = 0;
        for k = 1:nBlock
          R{k} = obj.getInterpolation(f, codim, {k});
          if nBlock>1
            fprintf(repmat('\b',1,length(s)));
            s = sprintf('progress evalInterpolation: %d / %d', k, nBlock); fprintf(s);
          end
        end
        if nBlock>1, fprintf('\n'); end
        R = cell2mat(R);
        val = sum(abs(R)>0,2);
        val(val==0) = 1;
        R = sum(R,2)./val;
      else
        basis = obj.evalGlobalBasis([], codim, 0, varargin{:}); % nExnBxnPxnC
        F = permute(obj.evalFunction(f, [], codim, [], varargin{:}), [1 4 2 3]); % nEx1xnPxnC
        [~,~,jac] = obj.evalTrafoInfo([], codim, varargin{:}); % nExnP
        [~, weights] = obj.getQuadData(codim); % nPx1
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
        lhs = sum(bsxfun(@times, permute(basis,[1 2 5 3 4]), ...
                                 permute(basis,[1 5 2 3 4])), 5); % nExnBxnBxnP
        lhs = permute(sum(bsxfun(@times, lhs, permute(dX, [1 3 4 2])), 4),[2 3 1]); % nBxnBxnE
        rhs = sum(bsxfun(@times, F, basis), 4); % nExnBxnP
        rhs = sum(bsxfun(@times, rhs, permute(dX, [1 3 2])),3).'; % nBxnE
        dMap = reshape(obj.getDoFMap(codim, varargin{:}),[],1); % nB*nE
        iZ = (dMap~=0); dMap = abs(dMap(iZ)); % 'nB*nE'
        lhs = blockify(lhs); % (nB*nE)x(nB*nE)
        rhs = accumarray(dMap,lhs(iZ,iZ)\rhs(iZ))./accumarray(dMap,1);
        I = unique(dMap);
        R = zeros(obj.getNDoF(),1); % nDoFx1
        R(I) = rhs(I);
      end
    end
  end
end
