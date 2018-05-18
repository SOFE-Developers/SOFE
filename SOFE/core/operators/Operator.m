classdef Operator < SOFE
  properties
    codim = 0;
    data, dataCache
    fesTrial, fesTest
    matrix
    matrix0
    loc, idx
    state
  end
  methods % constructor
    function obj = Operator(data, fesTrial, varargin) % [fesTest loc]
      if isnumeric(data) && numel(data)<4
        data = @(x)data+zeros(size(x,1),numel(data)); 
      end
      obj.dataCache = data;
      obj.data = data;
      obj.fesTrial = fesTrial;
      obj.fesTrial.register(obj);
      if nargin > 2 && ~isempty(varargin{1})
        obj.fesTest = varargin{1};
        obj.fesTest.register(obj);
        % sync quadRules
        if obj.fesTrial.quadRule{1}.order > obj.fesTest.quadRule{1}.order
          obj.fesTest.setQuadRule(obj.fesTrial.quadRule);
        else
          obj.fesTrial.setQuadRule(obj.fesTest.quadRule);
        end
      else
        obj.fesTest = fesTrial;
      end
      obj.idx = ':';
      if nargin > 3
        obj.loc = varargin{2};
      end
      try
        obj.matrix0 = obj.assembleRef();
      catch
      end
    end
    function notify(obj, varargin) % [time, state, dState]
      if nargin < 2
        obj.matrix = [];
        obj.idx = ':';
        obj.notifyObservers();
      else
        try, obj.state.U =  varargin{2}; end
        try, obj.state.dU =  varargin{3}; end
        if ~isnumeric(obj.dataCache)
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
        if ~isempty(obj.loc)
          if nargin(obj.loc) > 1 % loc(x,t)
            obj.matrix = [];
            obj.idx = obj.fesTrial.mesh.isBoundary(@(x)obj.loc(x, varargin{1}));
          else
            if strcmp(obj.idx, ':')
              obj.idx = obj.fesTrial.mesh.isBoundary(@(x)obj.loc(x));
            end
          end
        end
      end
    end
  end
  methods % assemble % apply
    function assemble(obj)
      % single sparse assembly
      if ~isempty(obj.matrix), return, end
      M = obj.fesTest.getNDoF(); N = obj.fesTrial.getNDoF();
      obj.matrix = sparse(M, N);
      if ~any(obj.idx), return, end
      nBlock = obj.fesTrial.nBlock(obj.codim+1);
      for k = 1:nBlock
        e = []; r = []; c = [];
        I = obj.fesTrial.getBlock(obj.codim, k);
        if ~isempty(I)
          e = obj.assembleOp(k); % nExnBxnB
          r = obj.fesTest.getDoFMap(obj.codim, I); % nBxnE
          c = obj.fesTrial.getDoFMap(obj.codim, I); % nBxnE
          r = repmat(abs(r)',[1 1 size(c,1)]); % nExnBxnB
          c = permute(repmat(abs(c)',[1 1 size(r,2)]), [1 3 2]); % nExnBxnB
          if ~ischar(obj.idx)
            e=e(obj.idx(I),:,:); r=r(obj.idx(I),:,:); c=c(obj.idx(I),:,:);
          end
        end
        I = (r.*c==0); if any(I(:)), r(I) = []; c(I) = []; e(I) = []; end %#ok<AGROW>
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
    end
    function R = integrate(obj, hasCoeff, basisI, basisJ, k)
      [~, weights] = obj.fesTrial.getQuadData(obj.codim);
      [~,~,jac] = obj.fesTrial.evalTrafoInfo([], obj.codim, {k}); % nExnP
      if hasCoeff
        if isnumeric(obj.data)
          coef = obj.fesTrial.evalDoFVector(obj.data, [], obj.codim, 0, {k}); % nExnPx(nD*nD)
        else
          try, S = obj.observers{1}.evalState(k); catch, S = obj.state; end
          coef = obj.fesTrial.evalFunction(obj.data, [], obj.codim, S, {k}); % nExnP
        end
      else
        coef = 1.0;
      end
      dX = bsxfun(@times, abs(jac).*coef, weights'); % nExnP
      nE = size(basisI, 1); nBI = size(basisI, 2); nBJ = size(basisJ,2); nP = size(basisI,3);
      try
        tprod(1,1,1,1);
        basisI = bsxfun(@times, basisI, permute(dX,[1 3 2]));
        R = tprod(reshape(basisI,nE,nBI,[]), ...
                  reshape(basisJ,nE,nBJ,[]), [1 2 -1], [1 3 -1]);
      catch err
        fprintf(['tprod:' err.message '\n']);
        basisI = reshape(basisI, nE, nBI, nP, []);
        basisJ = reshape(basisJ, nE, nBJ, nP, []);
        R = sum(bsxfun(@times, permute(basisI, [1 2 5 3 4]), ...
                               permute(basisJ, [1 5 2 3 4])), 5); % nExnBIxnBJxnPx(nC*nD)
        R = sum(bsxfun(@times, R, permute(dX, [1 3 4 2])), 4); % nExnBIxnBJ
      end
    end
    function R = apply(obj, x)
      R = obj.matrix*x; return
      % assemble on the fly
      dm = obj.fesTrial.getDoFMap(0); %#ok<UNRCH>
      R = accumarray(dm(:),reshape(obj.matrix0*x(dm),[],1));
    end
  end
end
