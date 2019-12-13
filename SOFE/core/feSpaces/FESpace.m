classdef FESpace < SOFE
% FESpace Class representing the Finite Element Space
%
%   constructor
%   -----------
%     obj = FESpace(mesh, element [, fixB, shift])
%
%     where
%     
%     mesh ... (sub)class of type Mesh
%     element --- (sub)class of type Element
%     fixB ... function handle for location of Dirichlet boundary
%     shift ... function handle for values of Dirichlet boundary
%
%   properties
%   ----------
%     mesh 
%     element
%     fixB
%     shift
%     freeDoFs
% 
%   public methods
%   --------------
%     evalFunction / evalReferenceMap / evalTrafoInfo
%     evalGlobalBasis / evalDoFVector
%     getDoFMap / getNDoF / extractDoFs / getBoundaryDoFs / getFreeDoFs
%     getL2Projection / getL2Interpolant / getLagrangeInterpolant
%     evalJumpResidual / getRecoveredGradient
%
%     See also meshes/Mesh, elements/Element
%
% SOFE Toolbox.
% Copyright 2018, Dr. Lars Ludwig
  properties
    mesh 
    element
    fixB
    shift
    narginShift
    freeDoFs
    cache
    isCaching = false;
    nBlock
    nBlockGlobal
  end
  methods % constructor, observer and caching.
    function obj = FESpace(mesh, element, varargin) % [fixB, shift]
      try
        if mesh.element.isSimplex ~= element.isSimplex || ...
         mesh.element.dimension ~= element.dimension
          warning('! Mesh and Element are not compatible !');
        end
      catch
      end
      obj.element = element;
      obj.mesh = mesh; obj.mesh.register(obj);
      obj.setBlocking();
      obj.fixB = @(x)false(size(x,1),1);
      obj.shift = []; obj.narginShift = 0;
      if nargin > 2
        obj.fixB = varargin{1};
        if nargin > 3 && ~isempty(varargin{2})
          obj.shift = varargin{2};
          try
            obj.narginShift = nargin(obj.shift);
          catch
            obj.narginShift = 0; % scalar shift
          end
        end
      end
    end
    function notify(obj, varargin) % [message]
      obj.setBlocking(); % resets cache
    end
    function setElement(obj, element)
      obj.element = element;
      obj.notify('elementChanged');
    end
    function resetCache(obj)
      nBlockMax = max(obj.nBlock);
      obj.cache.Phi = cell(nBlockMax, 1);
      obj.cache.DPhiInv = cell(nBlockMax, 1);
      obj.cache.jac = cell(nBlockMax, 1);
      obj.cache.basis = cell(nBlockMax, 1);
      for k = 1:nBlockMax
        obj.cache.Phi{k} = cell(obj.element.dimension+1,3);
        obj.cache.DPhiInv{k} = cell(obj.element.dimension+1,1);
        obj.cache.jac{k} = cell(obj.element.dimension+1,1);
        obj.cache.basis{k} = cell(obj.element.dimension+1,3);
      end
      obj.cache.dM = [];
      obj.cache.doFMap = cell(nBlockMax,1);
      for k = 1:nBlockMax
        obj.cache.doFMap{k} = cell(obj.element.dimension+1,1);
      end
      obj.cache.nDoF = [];
      obj.freeDoFs = [];
      %
      obj.notifyObservers('fesChanged');
    end
    function setBlockingGlobal(obj, N)
      if isscalar(N), N = N*ones(obj.element.dimension+1,1); end
      obj.nBlockGlobal = N;
      obj.setBlocking();
    end
    function setBlocking(obj)
      if ~isempty(obj.nBlockGlobal)
        obj.nBlock = obj.nBlockGlobal;
      else
        nC = obj.element.getNC(); nD = obj.element.dimension;
        nB = obj.element.nB(end); nQ = numel(obj.element.quadRule{1}.weights);
        obj.nBlock = ones(nD+1,1);
        if nD > 1
          elPerBlock = max(1,SOFE.getElementsPerBlock(nB, nQ, nC, nD));
          obj.nBlock = ones(nD+1,1);
          obj.nBlock(1) = ceil(obj.mesh.topology.getNumber(nD)/elPerBlock);
          %
          nB = obj.element.nB(end-1); nQ = numel(obj.element.quadRule{2}.weights);
          elPerBlock = max(1,SOFE.getElementsPerBlock(nB, nQ, nC, nD));
          obj.nBlock(2) = ceil(obj.mesh.topology.getNumber(nD-1)/elPerBlock);
          %
          if nD > 2
            nB = obj.element.nB(end-2); nQ = numel(obj.element.quadRule{3}.weights);
            if nB>0
              elPerBlock = max(1,SOFE.getElementsPerBlock(nB, nQ, nC, nD));
              obj.nBlock(3) = ceil(obj.mesh.topology.getNumber(nD-2)/elPerBlock);
            end
          end
        end
      end
      obj.nBlock(end) = 1;
      obj.resetCache();
    end
    function R = getBlock(obj, codim, varargin) % [k]
      nE = obj.mesh.topology.getNumber(num2str(codim));
      if obj.nBlock(codim+1)>nE, error('!Number of blocks exceeds number of elements!'); end
      R = unique(floor(linspace(0,nE,obj.nBlock(codim+1)+1)));
      R = [R(1:end-1)+1; R(2:end)];
      if nargin > 2
        if varargin{1}>obj.nBlock(codim+1)
          R = [];
        else
          R = (R(1,varargin{:}):R(2,varargin{:}))';
        end
      else
        R = (1:R(end))';
      end
    end
  end
  methods % evaluation.
    function R = evalFunction(obj, F, points, codim, S, varargin) % [{k} or I]
      % R = evalFunction(obj, F, points, codim, S [,idx])
      %   evaluates function in local/quadrature points
      %
      % Input
      % -----
      % F = F(x [, u, gradU])
      %   function handle to be evaluated
      % points
      %   local evaluation points ([] if not needed)
      % codim
      %   codimension of quadrature points (in case of POINTS = [])
      % S.U
      %   nExnPxnC state array for nonlinear function F
      % S.dU
      %   nExnPxnCxnD gradient state array for nonlinear function F
      % idx
      %   subentities I or block {k}
      %
      % Output
      % ------
      % R
      %   nExnPxnC array
      %
      %--------------------------------
      % SOFE Toolbox.
      % Copyright 2017, Dr. Lars Ludwig
      %--------------------------------
      if isempty(points)
        points = obj.element.getQuadData(codim);
      else
        codim = obj.element.dimension-size(points,2);
      end
      points = {obj.evalReferenceMap(points, codim, 0, varargin{:})};
      R = obj.mesh.evalFunction(F, points, S, varargin{:});
    end
    function R = evalReferenceMap(obj, points, codim, order, varargin) % [{k} or I]
      % R = evalReferenceMap(obj, points, codim, order [, idx])
      %   evaluates reference map in local/quadrature points
      %
      % Input
      % -----
      % points
      %   local evaluation points ([] if not needed)
      % codim
      %   codimension of quadrature points (in case of POINTS = [])
      % idx
      %   subentities I or block {k}
      %
      % Output
      % ------
      % R
      %   nExnPxnW[xnD] array
      %
      %--------------------------------
      % SOFE Toolbox.
      % Copyright 2017, Dr. Lars Ludwig
      %--------------------------------
      if isempty(points)
        forcedPoints = obj.element.getQuadData(codim);
      else
        forcedPoints = points;
        codim = obj.element.dimension-size(points,2);
      end
      doCache = 0;
      if ~isempty(varargin) && iscell(varargin{1})
        doCache =  obj.isCaching && isempty(points);
        k = varargin{1};
        varargin{1} = obj.getBlock(codim, k{1});
      end
      %
      if doCache && ~isempty(obj.cache.Phi{k{1}}{codim+1,order+1})
        R = obj.cache.Phi{k{1}}{codim+1,order+1};
        return
      end
      R = obj.mesh.evalReferenceMap(forcedPoints, order, varargin{:});
      if doCache
        obj.cache.Phi{k{1}}{codim+1,order+1} = R;
      end
    end
    function R = evalNormalVector(obj, nFlag, varargin) % [{k} or I]
      if ~isempty(varargin) && iscell(varargin{1})
        varargin{1} = obj.getBlock(1, varargin{1}{1});
      end
      points = obj.element.getQuadData(1);
      R = obj.mesh.evalNormalVector(points, nFlag, varargin{:});
    end
    function [R, invR, jacR] = evalTrafoInfo(obj, points, codim, varargin) % [{k} or I]
      % [R, invR, jacR] = evalTrafoInfo(obj, points, codim [, idx])
      %   evaluates data for integral trafo in local/quadrature points
      %
      % Input
      % -----
      % points
      %   local evaluation points ([] if not needed)
      % codim
      %   codimension of quadrature points (in case of POINTS = [])
      % idx
      %   subentities I or block {k}
      %
      % Output
      % ------
      % R
      %   derivative of reference map, nExnPxnWxnD array
      % invR
      %   (pseudo) inverse of derivative of reference map, nExnPxnDxnW array
      % jacR
      %   (generalized) jacobian of reference map, nExnP array
      %
      %--------------------------------
      % SOFE Toolbox.
      % Copyright 2017, Dr. Lars Ludwig
      %--------------------------------
      if isempty(points)
        forcedPoints = obj.element.getQuadData(codim);
      else
        forcedPoints = points;
        codim = obj.element.dimension-size(points,2);
      end
      doCache = 0;
      if ~isempty(varargin) && iscell(varargin{1})
        doCache =  obj.isCaching && isempty(points);
        k = varargin{1};
        varargin{1} = obj.getBlock(codim, k{1});
      end
      %
      if doCache && ~isempty(obj.cache.Phi{k{1}}{codim+1,2})
        R = obj.cache.Phi{k{1}}{codim+1,1+1};
        invR = obj.cache.DPhiInv{k{1}}{codim+1};
        jacR = obj.cache.jac{k{1}}{codim+1};
        return
      end
      [R, invR, jacR] = obj.mesh.evalTrafoInfo(forcedPoints, varargin{:});
      if doCache
        obj.cache.Phi{k{1}}{codim+1,2} = R;
        obj.cache.DPhiInv{k{1}}{codim+1} = invR;
        obj.cache.jac{k{1}}{codim+1} = jacR;
      end
    end
    function R = computeGlobalBasis(obj, points, codim, order, varargin) % [{k} or I]
      if iscell(points)
        varargin{1} = points{2}; codim = obj.element.dimension - size(points{1},2);
        basis = obj.element.evalBasis(points{1}, order); % nBxnPxnC[xnD]
        pVec{2} = [2 3 1]; pVec{1} = [4 5 1 2 3];
      else
        if isempty(points)
          forcedPoints = obj.element.getQuadData(codim);
        else
          forcedPoints = points;
          codim = obj.element.dimension-size(points,2);
        end
        basis = obj.element.evalBasis(forcedPoints, order); % nBxnPxnC[xnD]
        pVec{2} = [1 3 2]; pVec{1} = [1 5 2 3 4];
      end
      if ~(order==0 && (strcmp(obj.element.conformity, 'H1') || strcmp(obj.element.conformity, 'L2')))
        [trafo{1},trafo{2},trafo{3}] = obj.evalTrafoInfo(points, codim, varargin{:}); % nExnP[xnWxnD]
        trafo{1} = permute(trafo{1}, pVec{1}); % nEx1xnPxnWxnD
        trafo{2} = permute(trafo{2}, pVec{1}); % nEx1xnPxnDxnW
        trafo{3} = permute(trafo{3}, pVec{2}); % nEx1xnP
      end
      basis = permute(basis, [6 1 2 3 4 5]); % 1xnBxnPxnC[xnD]
      switch obj.element.conformity
        case {'L2', 'H1'}
          switch order
            case 0
              R = basis; % 1xnBxnPxnC
            case 1
              try
                R = tprod(permute(basis,[2 3 4 5 1]), permute(trafo{2},[1 3 4 5 2]), [2 3 4 -1], [1 3 -1 5]); % nExnBxnPxnCxnD
              catch
                R = sum(bsxfun(@times, permute(basis,    [1 2 3 4 6 5]), ...
                                       permute(trafo{2}, [1 2 3 6 5 4])), 6); % nExnBxnPxnCxnD
              end
            case 2
              R = permute(basis, [1 2 3 6 4 5]); % % nExnBxnPxnD2xnCxnD1
              R = sum(bsxfun(@times, permute(R,    [1 2 3 4 5 7 6]), ...
                                     permute(trafo{2}, [1 2 3 7 6 5 4])), 7); % nExnBxnPxnD2xnCxnD1
              R = permute(R, [1 2 3 6 5 4]); % % nExnBxnPxnD1xnCxnD2
              R = sum(bsxfun(@times, permute(R,    [1 2 3 4 5 7 6]), ...
                                     permute(trafo{2}, [1 2 3 7 6 5 4])), 7); % nExnBxnPxnD2xnCxnD1
              R = permute(R, [1 2 3 5 4 6]); % % nExnBxnPxxnCxnD1xnD2
          end
        case 'HDiv'
          if codim == 0
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
          else % codim == 1
            N = obj.evalNormalVector(2, varargin{:}); % nExnPxnD
            switch order
              case 0
                R = basis;
              case 1
                R = sum(bsxfun(@times, permute(basis,    [1 2 3 4 6 5]), ...
                                       permute(trafo{2}, [1 2 3 6 5 4])), 6); % nExnBxnPx(nC=1)xnD
            end
            R = bsxfun(@times, R, permute(N, [1 4 2 3])); % nExnBxnPxnCx[...]
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
      S = 1-2*(dMap<0); % S = sign(dMap);
      R = bsxfun(@times, S, R); % nExnBxnPxnC[xnD]
    end
    function R = evalGlobalBasis(obj, points, codim, order, varargin) % [{k} or I]
      doCache = 0;
      if ~isempty(varargin) && iscell(varargin{1})
        doCache =  obj.isCaching && isempty(points);
        k = varargin{1};
      end      
      if doCache && ~isempty(obj.cache.basis{k{1}}{codim+1, order+1})
        R = obj.cache.basis{k{1}}{codim+1, order+1};
        return
      end
      if ~isempty(points)
        codim = obj.element.dimension-size(points,2);
      end
      if isempty(varargin) || (ischar(varargin{1}) && strcmp(varargin{1},':'))
        nBl = obj.nBlock(codim+1); R = cell(nBl,1);
        for k = 1:nBl
          R{k} = obj.evalGlobalBasis(points, codim, order, {k});
        end
        try
          R = cell2mat(R);
        catch
          R = padcell2mat(R);
        end
        [~,I] = sort(obj.getBlock(codim));
        R = R(I,:,:,:,:);
      else
        R = obj.computeGlobalBasis(points, codim, order, varargin{:});
      end
      if doCache
        obj.cache.basis{k{1}}{codim+1, order+1} = R;
      end
    end
    function R = evalDoFVector(obj, U, points, codim, order, varargin) % [{k} or I]
      assert(size(U,1)==obj.getNDoF(), 'First argument must be DoFVector(s)!');
      if iscell(points)
        R = obj.evalDoFVectorGlobal(U, points, order);
      else
        R = obj.evalDoFVectorLocal(U, points, codim, order, varargin{:});
      end
    end
    function R = evalDoFVectorGlobal(obj, U, points, order)
      if numel(points) == 1
        points = obj.mesh.evalInversReferenceMap(points{1});
      end
      isValid = points{2}>0;
      if ~any(isValid), R = nan(size(isValid)); return; end
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
      if isempty(varargin), varargin{1} = ':'; end
      if ~isempty(points)
        codim = obj.element.dimension-size(points,2);
      end
      if iscell(varargin{1})
        assert(codim==0 || order==0, '! Derivatives for traces not supported !');
        basis = obj.evalGlobalBasis(points, codim, order, varargin{:}); % [1/nE]xnB[xnP]xnCx[nD]x[nD]
        dMap = abs(obj.getDoFMap(codim, obj.getBlock(codim, varargin{1}{1}))).'; % nExnB
        if isempty(dMap), R = []; return; end
        R = zeros(numel(dMap), size(U,2)); I = dMap > 0; R(I,:) = U(dMap(I),:); % (nE*nB)xN
        R = reshape(R, size(dMap,1), size(dMap,2), []); % nExnBxN
        R = sum(bsxfun(@times,permute(R,[1 4:7 3 2]),permute(basis,[1 3:7 2])),7); % nExnPxnCx[nD]x[nD]xN
        sz = [size(R) ones(1,6-ndims(R))]; 
        R = reshape(permute(R, [1 2 3 6 4 5]), [sz([1 2])  sz(3)*sz(6) sz([4 5])]); % nExnPx[nC*N]x[nD]x[nD]
      else % all blocks
        nBl = obj.nBlock(codim+1); R = cell(nBl,1); s = '';
        for k = 1:nBl
          R{k} = obj.evalDoFVectorLocal(U, points, codim, order, {k}); % nExnPxnCxnD
          if nBl>1 
            fprintf(repmat('\b',1,length(s)));
            s = sprintf('progress evalDoFVector: %d / %d', k, nBl); fprintf(s);
          end
        end
        if nBl>1, fprintf('\n'); end
        try R = cell2mat(R); catch, R = padcell2mat(R); end
        R = R(varargin{1},:,:,:,:);
      end
    end
  end
  methods % dof mapping
    function R = getNDoF(obj)
      if isempty(obj.cache.nDoF)
        obj.cache.nDoF = sum(prod(obj.element.doFTuple,1)'.*obj.mesh.topology.getNumber());
      end
      R = obj.cache.nDoF;
    end
    function R = getDoFMap(obj, codim, varargin) % [{k} or I]
      if isempty(varargin), varargin{1} = ':'; end
      if iscell(varargin{1})
        if ~isempty(obj.cache.doFMap{varargin{1}{1}}{codim+1})
          R = obj.cache.doFMap{varargin{1}{1}}{codim+1};
          return
        end
        k = varargin{1}; varargin{1} = obj.getBlock(codim, k{1});
        dim = obj.mesh.topology.dimP - codim;
        doFTuple = prod(obj.element.doFTuple,1)';
        csnDoF = [0; cumsum(doFTuple.*obj.mesh.topology.getNumber())];
        R = cell(dim+1,1); % {nD+1}((nBLoc*nESub)xnE)
        for d = 0:dim
          C = obj.mesh.topology.connectivity{dim+1, d+1}; % nE(dim)xnE(d)
          if doFTuple(d+1) ~= 0
            R{d+1} = csnDoF(d+1) + (C(varargin{1},:)'-1)*doFTuple(d+1);
            R{d+1} = kron(R{d+1},ones(doFTuple(d+1),1)) + kron(ones(size(R{d+1})),(1:doFTuple(d+1))');
            R{d+1} = obj.orient(R{d+1}, dim, d, varargin{:}); % {nD+1}((nBLoc*nESub)*nE)
          end
        end
        R = cell2mat(R); % nBxnE
        obj.cache.doFMap{k{1}}{codim+1} = R;
      else
        nBl = obj.nBlock(codim+1);
        R = cell(nBl,1);
        for k = 1:nBl
          R{k} = obj.getDoFMap(codim, {k}); % nExnPxnCxnD
        end
        R = cell2mat(R');
        if ~isempty(R) && isreal(varargin{1})
          R = R(:,varargin{1});
        end
      end
    end
    function R = orient(obj, D, dim, d, varargin) % [I]
      orient = permute(obj.mesh.topology.getOrientation(dim, d, varargin{:}),[2 1 3]); % nESubxnExnO
      if isempty(varargin), nE = obj.mesh.topology.getNumber(dim); else, nE = numel(varargin{1}); end
      if isempty(orient)
        R=D;
      else 
        nEntSub = obj.element.getNEntSub(dim);
        D = reshape(D, [], nEntSub(d+1)*nE); % nBLocx(nESub*nE)
        orient = reshape(orient, [], size(orient, 3)); % (nESub*nE)xnO
        if isa(obj.element, 'PpH1') && dim==3 && d==2
          R = zeros(3*size(D,1),size(D,2)); % nBLocx(nESub*nE)
        else
          R = zeros(size(D)); % nBLocx(nESub*nE)
        end
        for k = 1:obj.mesh.topology.nO
          iO = orient==k; % (nESub*nE)
          if any(iO)
            doFEnum = obj.element.getDoFEnum(d,k); % nBLoc
            if ~isempty(doFEnum)
              R(abs(doFEnum),iO) = bsxfun(@times, sign(doFEnum(:)), D(:, iO)); % nBLocx(nESub*nE)
            end
          end
        end
      end
      %
      if dim == obj.element.dimension && d == dim-1 && strcmp(obj.element.conformity, 'HDiv')
        orient = reshape(obj.mesh.topology.getNormalOrientation(varargin{:}).',1,[]); % 1x(nESub*nE)
        R = bsxfun(@times, R, orient); % nBLocx(nESub*nE)
      end
      %
      R = reshape(R, [], nE); % (nBLoc*nESub)xnE
    end
    function R = extractDoFs(obj, codim, I)
      dMap = obj.getDoFMap(codim);
      nC = size(I,2);
      doFs = cell(nC,1);
      for d = 1:nC
        doFs{d} = reshape(abs(dMap(d:nC:end, I(:,d))), [], 1);
      end
      R = accumarray(setdiff(unique(reshape(cell2mat(doFs),[],1)), 0), 1, [obj.getNDoF() 1])>0;
    end
    function R = getBoundaryDoFs(obj, varargin) % [loc]
      R = obj.extractDoFs(1, obj.mesh.isBoundary(varargin{:}));
    end
    function R = getFreeDoFs(obj)
      if ~isempty(obj.freeDoFs)
        R = obj.freeDoFs;
      else
        if obj.element.dimension > 0
          R = ~obj.getBoundaryDoFs(obj.fixB);
        else
          R = true(obj.getNDoF(),1);
        end
        obj.freeDoFs = R;
      end
    end
    function R = getProjector_(obj)
      R = obj.mesh.topology.getProjector();
      if ~isempty(R), assert(obj.element.order==1, 'TODO: higher order projector'); end
    end
  end
  methods % refinement
    function R = uniformRefine(obj)
      if nargout>0
        dM0 = obj.getDoFMap(0);
      end
      obj.mesh.uniformRefine();
      if nargout>0
        R = obj.getProlongator(dM0, obj.getDoFMap(0));
      end
    end
  end
  methods % interpolation
    function R = getShift(obj, varargin) % [time]
      R = zeros(obj.getNDoF,1);
      if ~isempty(obj.shift)
        if obj.narginShift > 1
          if isempty(varargin)
            func = @(x)obj.shift(x, 0.0);
          else
            func = @(x)obj.shift(x, varargin{1});
          end
        else
          func = obj.shift;
        end
        R = obj.getL2Interpolant(func, obj.mesh.element.dimension-1, obj.mesh.isBoundary()); % nDoFx1
      end
    end
    function R = getL2Projection(obj, f)
      mass = OpIdId(1, 0, obj); mass.assemble();
      l2 = FcId(f, obj, 0); l2.assemble();
      R = mass.matrix\l2.matrix;
    end
    function R = getL2Interpolant(obj, f, dim, varargin) % [I]
      codim = obj.element.dimension - dim;
      if isempty(varargin), varargin{1} = ':'; end
      if dim==0
        R = zeros(obj.getNDoF(),1);
        dMap = obj.getDoFMap(codim, varargin{1}); % nBxnE
        if isempty(dMap), return; end
        R(dMap') = f(obj.mesh.nodes(varargin{1},:));
      else
        I = unique(obj.mesh.topology.connectivity{dim+1, dim}(varargin{1},:));
        R = obj.getL2Interpolant(f, dim-1, I); % recursion
        %
        nEntSub = obj.element.getNEntSub(dim);
        doFTuple = prod(obj.element.doFTuple(:,1:dim),1);
        if dim==3 && obj.element.isSimplex && ~obj.element.isLagrange
          doFTuple(3) = 3*doFTuple(3);
        end
        if dim==3 && isa(obj.element, 'NdPp')
          doFTuple(3) = 1.5*doFTuple(3);
        end
        if dim==3 && isa(obj.element, 'BDMPRot')
          doFTuple(3) = 1.5*(doFTuple(3)-3*(doFTuple(2)-2)) + 3*(doFTuple(2)-2);
        end
        offsetB = sum(doFTuple.*nEntSub(1:dim));
        dMap = obj.getDoFMap(codim, varargin{:}); % nBxnE
        dMap = reshape(dMap(offsetB+1:end,:),[],1); % nB*nE
        if isempty(dMap), return; end
        inZ = (dMap~=0); dMap = abs(dMap(inZ)); % 'nB*nE
        %
        F = obj.evalFunction(f,[],codim,[],varargin{:}); % nExnPxnC
        F = bsxfun(@minus,F,obj.evalDoFVector(R,[],codim,0,varargin{:}));
        [~,~,jac] = obj.evalTrafoInfo([], codim, varargin{:}); % nExnP
        [~, weights] = obj.element.getQuadData(codim); % nPx1
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
        basis = obj.evalGlobalBasis([], codim, 0, varargin{:}); % nExnBxnPxnC
        basis = basis(:,offsetB+1:end,:,:); % nExnBxnPxnC
        lhs = sum(bsxfun(@times, permute(basis,[1 2 5 3 4]), ...
                                 permute(basis,[1 5 2 3 4])), 5); % nExnBxnBxnP
        lhs = permute(sum(bsxfun(@times, lhs, permute(dX, [1 3 4 2])), 4),[2 3 1]); % nBxnBxnE
        rhs = sum(bsxfun(@times, permute(F, [1 4 2 3]), basis), 4); % nExnBxnP
        rhs = sum(bsxfun(@times, rhs, permute(dX, [1 3 2])),3).'; % nBxnE
        lhs = blockify(lhs); % (nB*nE)x(nB*nE)
        %
        R(dMap) = lhs(inZ,inZ)\reshape(rhs(inZ),[],1);
      end
    end
    function R = getLagrangeInterpolant(obj, f, dim, varargin) % [I]
      assert(obj.element.isLagrange, 'Element basis must have Lagrange property');
      codim = obj.element.dimension - dim;
      if isempty(varargin), varargin{1} = ':'; end
      if dim==0
        R = zeros(obj.getNDoF(),1); % nDoFx1
        dMap = obj.getDoFMap(codim, varargin{1}); % nBxnE
        if isempty(dMap), return; end
        R(dMap') = f(obj.mesh.nodes(varargin{1},:)); % nDoFx1
      else
        I = unique(obj.mesh.topology.connectivity{dim+1, dim}(varargin{1},:));
        R = obj.getLagrangeInterpolant(f, dim-1, I); % recursion
        %
        nEntSub = obj.element.getNEntSub(dim);
        doFTuple = prod(obj.element.doFTuple(:,1:dim),1);
        if dim==3 && isa(obj.element, 'NdPp')
          doFTuple(3) = 1.5*doFTuple(3);
        end
        if dim==3 && isa(obj.element, 'BDMPRot')
          doFTuple(3) = 1.5*(doFTuple(3)-3*(doFTuple(2)-2)) + 3*(doFTuple(2)-2);
        end
        offsetB = sum(doFTuple.*nEntSub(1:dim));
        dMap = obj.getDoFMap(codim, varargin{:}); % nBxnE
        dMap = reshape(dMap(offsetB+1:end,:),[],1); % nB*nE
        if isempty(dMap), return; end
        [points, dir] = obj.element.getLagrangePoints(dim, obj.element.order); % nPxnD
        points = points(sum(prod(obj.element.doFTuple(:,1:dim),1).*nEntSub(1:dim))+1:end,:); % nPxnD
        dir = dir(sum(prod(obj.element.doFTuple(:,1:dim),1).*nEntSub(1:dim))+1:end,:); % nPxnD
        P = obj.mesh.evalReferenceMap(points, 0, varargin{1}); % nExnPxnD
        F = f(reshape(P, [], size(P,3))); % (nE*nP)xnC
        switch obj.element.conformity
          case {'L2', 'H1'}
            D = repmat(permute(eye(obj.element.getNC),[3 4 1 2]),size(P,1),size(P,2));
          case 'HDiv'
            switch codim
              case 1
                N = reshape(obj.mesh.evalNormalVector(points, 0, varargin{1}),[],size(F,2)); % (nE*nP)xnD
                R(dMap(:)) = reshape(dot(F,N,2), size(P,1),[])'; % nDoFx1
                return
              case 0
                [~, D, J] = obj.mesh.evalTrafoInfo(points); % nExnPxnDxnD
                D = J.*permute(D, [1 2 4 3]); % nExnPxnDxnD
            end
          case 'HRot'
            D = obj.mesh.evalReferenceMap(points, 1, varargin{1}); % nExnPxnCxnD
        end
        D = sum(bsxfun(@times, D, permute(dir, [3 1 4 2])),4);
        R(dMap(:)) = reshape(dot(F,reshape(D,size(F)), 2), size(P,1), [])'; % nDoFx1
      end
    end
  end
  methods % DoFVector operations
    function R = evalJumpResidual(obj, U, points, order)
      assert(obj.mesh.dimW==2, '!Dim of world equals 2 assumed!');
      %
      nOrient = obj.mesh.topology.getNormalOrientation();
      e2F = obj.mesh.topology.getElem2Face(); % nEx3
      if order<1, nD = 1; else, nD = obj.mesh.topology.dimP; end
      R = zeros(obj.mesh.topology.getNumber('1'), 2, numel(points), obj.element.getNC(), nD); % nFx(L/R)xnPxnCxnD
      for fLoc = 1:3
        for k = 1:2
          P = obj.mesh.topology.upliftPointsN(reshape(points,[],1), fLoc, k);
          tmp = obj.evalDoFVector(U, P, [], order); % nExnPxnCxnD
          I = nOrient(:,fLoc)==(3-2*k);
          R(e2F(I,fLoc),k,:,:,:) = tmp(I,:,:,:);
        end
      end
      R = permute(R, [1 3 4 5 2]); % nFxnPxnCxnDx(L/R)
    end
    function R = getRecoveredGradient(obj, U)
      nD = obj.element.dimension; nB = obj.element.nB(nD);
      points = linspace(0,1,obj.element.order+1)'; P = points;
      for i = 2:nD
        P = [kron(ones(length(points),1),P), kron(points,ones(length(points)^(i-1),1))];
      end
      if obj.element.isSimplex
        P = P((sum(P,2)<=1),:);
      end
      lhs = reshape(obj.element.evalBasis(P,0), nB, []).'; % (nP*nC)xnB
      dU = obj.evalDoFVector(U, P, [], 1); % nExnPxnCxnD
      dMap = obj.getDoFMap(0); % nBxnE
      R = zeros(obj.getNDoF,nD);
      for d = 1:nD
        rhs = sign(dMap).*(lhs\reshape(permute(dU(:,:,:,d), [2 3 1]), nB, [])); % nBxnE
        R(:,d) = accumarray(dMap(:),rhs(:))./accumarray(dMap(:),1); 
      end
    end
  end
  methods % prolong & restrict
    function R = getProlongator(obj, dM0, dM1, varargin) % [iE, iCh]
      try iE = varargin{1}; catch, iE = true(size(dM0,2),1); end
      try iCh = varargin{2}; catch, iCh = reshape(1:2^obj.element.dimension*numel(iE),numel(iE),[]); end
      nCh = size(iCh,2); nB = size(dM0,1);
      pCoeff = obj.element.getProlongationCoeff();
      I = cell(nCh,1); J = cell(nCh,1); C = cell(nCh,1);
      for k = 1:nCh
        I{k} = kron(ones(1,nB),dM1(:, iCh(:,k))'); % child
        J{k} = kron(dM0(:,iE)',ones(1,nB)); % parent
        C{k} = repmat(reshape(pCoeff(:,:,k),1,[]), numel(iE), 1);
      end
      I = cell2mat(I); J = cell2mat(J); C = cell2mat(C);
      C = C.*sign(I).*sign(J); I = abs(I); J = abs(J);
%       I = I(I~=0); J = J(I~=0); C = C(I~=0);
      cnt = sparse(I,J,1);
      R = sparse(I,J,C);
      [ii,jj,cnt] = find(cnt);
      idx = ii+size(R,1)*(jj-1);
      R(idx) = R(idx)./cnt;
%       R = bsxfun(@rdivide, sparse(I,J,C), accumarray(I(:),1)/nB); % for non-conforming elements
    end
    function R = prolongDoFVector(obj, U,dM0, dM1)
      nDoF1 = max(abs(dM1(:)));
      pCoeff = obj.element.getProlongationCoeff();
      R = prolongIdx(pCoeff, U, dM0, dM1, nDoF1); % --> TODO: prolongDoFVector.mex
    end
    function R = restrictDoFVector(obj, U, dM0, dM1)
      nDoF0 = max(abs(dM0(:)));
      nDoF1 = max(abs(dM1(:)));
      degDoF1 = accumarray(abs(dM1(:)), 1, [nDoF1, 1]);
      pCoeff = obj.element.getProlongationCoeff();
      R = restrictIdx(pCoeff, U, dM0, dM1, degDoF1, nDoF0); % --> restrictDoFVector.mex (degDoF incr in loop)
    end
  end
end