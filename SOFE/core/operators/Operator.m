classdef Operator < SOFE
  properties
    codim
    data, dataCache
    fesTrial, fesTest
    matrix
    loc
    hasCoeff = true
    matrixFree = false
    A0
  end
  methods % constructor
    function obj = Operator(data, fesTrial, varargin) % [fesTest loc]
      if isnumeric(data) && size(data,1)==1
        obj.dataCache = @(x)data+zeros(size(x,1),numel(data)); 
      else
        obj.dataCache = data;
      end
      obj.codim = 0;
      obj.data = data;
      obj.fesTrial = fesTrial;
      obj.fesTrial.register(obj);
      if ~isempty(varargin) && ~isempty(varargin{1})
        obj.fesTest = varargin{1};
        obj.fesTest.register(obj);
        obj.syncQuadRules();
        nBlock = max(obj.fesTrial.nBlock, obj.fesTest.nBlock);
        obj.fesTrial.setBlockingGlobal(nBlock);
        obj.fesTest.setBlockingGlobal(nBlock);
      else
        obj.fesTest = fesTrial;
      end
      if nargin > 3
        obj.loc = varargin{2};
      end
    end
    function syncQuadRules(obj)
      if obj.fesTrial.element.quadRule{1}.order > obj.fesTest.element.quadRule{1}.order
        obj.fesTest.element.quadRule = obj.fesTrial.element.quadRule;
      else
        obj.fesTrial.element.quadRule = obj.fesTest.element.quadRule;
      end
    end
    function notify(obj, message, varargin) % [time]
      switch message
        case 'fesChanged'
          if ~obj.matrixFree
            obj.matrix = [];
          end
          obj.notifyObservers('opChanged');
        case 'stateChanged'
          if ~isnumeric(obj.dataCache) % isFunctionHandle
            if nargin(obj.dataCache) == 2 % f(x,t)
              obj.matrix = [];
              obj.data = @(x)obj.dataCache(x, varargin{1});
            elseif nargin(obj.dataCache) == 3 % f(x,t,U)
              obj.matrix = [];
              obj.data = @(x, U)obj.dataCache(x, varargin{1}, U);
            elseif nargin(obj.dataCache) == 4 % f(x,t,U,d)
                obj.matrix = [];
                obj.data = @(x, U, d)obj.dataCache(x, varargin{1}, U, d);
            end
          end
        otherwise
          error('Unknown message');
      end
    end
  end
  methods % assemble, integrate & apply
    function assemble(obj)
      if ~isempty(obj.matrix), return, end
      M = obj.fesTest.getNDoF(); N = obj.fesTrial.getNDoF();
      obj.matrix = sparse(M, N);
      if ~isempty(obj.loc)
        idx = obj.fesTrial.mesh.isBoundary(@(x)obj.loc(x));
      else
        idx = ':';
      end
      if ~any(idx), return, end
      nBlock = obj.fesTrial.nBlock(obj.codim+1);
      obj.A0 = cell(nBlock,1);
      for k = 1:nBlock
        e = []; r = []; c = [];
        I = obj.fesTrial.getBlock(obj.codim, k);
        if ~isempty(I)
          e = obj.assembleOp(k); % nExnBxnB
          r = obj.fesTest.getDoFMap(obj.codim, {k}); % nBxnE
          c = obj.fesTrial.getDoFMap(obj.codim, {k}); % nBxnE
          r = repmat(abs(r)',[1 1 size(c,1)]); % nExnBxnB
          c = permute(repmat(abs(c)',[1 1 size(r,2)]), [1 3 2]); % nExnBxnB
          if ~strcmp(idx,':')
            e = e(idx(I),:,:);
            r = r(idx(I),:,:);
            c = c(idx(I),:,:);
          end
        end
        if obj.matrixFree
          sgnTest = sign(obj.fesTest.getDoFMap(obj.codim, {k}))'; % nExnB
          sgnTrial = permute(sign(obj.fesTrial.getDoFMap(obj.codim, {k}))', [1 3 2]); % nExnB
          obj.A0{k} = (e.*sgnTrial).*sgnTest; % nBxnBxnE
        end
        I = (r.*c==0); if any(I(:)), r(I) = []; c(I) = []; e(I) = []; end %#ok<AGROW>
        %
        try
          fsparse(0);
          obj.matrix = obj.matrix + fsparse(r(:), c(:), e(:), [M, N]);
        catch
          obj.matrix = obj.matrix + sparse(r(:), c(:), e(:), M, N);
        end
        if k>1
          if k>2
            fprintf(repmat('\b',1,length(s)));
          end
          s = sprintf(['progress assembly ', class(obj), ': %d / %d'], k, nBlock);
          fprintf(s);
        end
      end
      if k>1, fprintf('\n'); end
      obj.A0 = permute(cell2mat(obj.A0), [2 3 1]); % nBxnBxnE
    end
    function R = integrate(obj, basisI, basisJ, k)
      [~, weights] = obj.fesTrial.element.getQuadData(obj.codim);
      [~,~,jac] = obj.fesTrial.evalTrafoInfo([], obj.codim, {k}); % nExnP
      if obj.hasCoeff
        if isnumeric(obj.data)
          if size(obj.data,1)==1
            coef = permute(obj.data, [3 1 2]); % 1x1xnC
          else
            coef = obj.fesTrial.evalDoFVector(obj.data, [], obj.codim, 0, {k}); % nExnPxnC
          end
        else
          try S = obj.observers{1}.evalState(k); catch, S = []; end
          coef = obj.fesTrial.evalFunction(obj.data, [], obj.codim, S, {k}); % nExnPxnC
        end
      else
        coef = 1.0;
      end
      dX = bsxfun(@times, abs(jac).*coef, weights'); % nExnP
      nE = size(basisI, 1); nBI = size(basisI, 2); nBJ = size(basisJ,2); nP = size(basisI,3);
      try
        tprod(1,1,1,1);
        basisI = bsxfun(@times, basisI, permute(dX,[1 3 2])); % nExnBIxnPxnCxnD
        R = tprod(reshape(basisI,nE,nBI,[]), ...
                  reshape(basisJ,nE,nBJ,[]), [1 2 -1], [1 3 -1]);
      catch err
        fprintf(['tprod:' err.message '\n']);
        basisI = reshape(basisI, nE, nBI, nP, []); % nExnBIxnPx(nC*nD)
        basisJ = reshape(basisJ, nE, nBJ, nP, []); % nExnBJxnPx(nC*nD)
        R = sum(bsxfun(@times, permute(basisI, [1 2 5 3 4]), ...
                               permute(basisJ, [1 5 2 3 4])), 5); % nExnBIxnBJxnP
        R = sum(bsxfun(@times, R, permute(dX, [1 3 4 2])), 4); % nExnBIxnBJ
      end
    end
    function R = apply(obj, x, varargin) % [freeI, freeJ]
      if ~isempty(varargin)
        X = zeros(obj.fesTrial.getNDoF(), 1);
        X(varargin{2}) = x;
      else
        X = x;
      end
      dM = obj.fesTrial.getDoFMap(0);
      R = obj.A0*X(dM);
      dM = obj.fesTest.getDoFMap(0);
      R = accumarray(dM(:),R(:));
      if ~isempty(varargin)
        R = R(varargin{1});
      end
    end
    function R = apply_(obj, x, varargin) % [freeI, freeJ] % deprecated
      freeI = ':'; freeJ = ':';
      if ~isempty(varargin), freeI = varargin{1}; freeJ = varargin{2}; end
      X = zeros(obj.fesTrial.getNDoF(), 1);
      X(freeJ) = x;
      R = obj.matrix*X;
      R = R(freeI);
    end
  end
end
