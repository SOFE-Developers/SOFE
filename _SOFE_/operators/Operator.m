classdef Operator < SOFEClass
  properties
    codim = 0;
    data, dataCache
    fesTrial, fesTest
    matrix
    state
    loc, idx
  end
  methods % constructor
    function obj = Operator(data, fesTrial, varargin) % [fesTest loc]
      if isreal(data), data = @(x)data+zeros(size(x,1),numel(data)); end
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
    end
    function notify(obj, varargin) % [time, state]
      if nargin < 2
        obj.matrix = [];
        obj.idx = ':';
      else
        if nargin(obj.dataCache) == 2 % f(x,t)
          obj.matrix = [];
          obj.data = @(x)obj.dataCache(x, varargin{1});
        elseif nargin(obj.dataCache) == 3 % f(x,t,U)
          obj.matrix = [];
          obj.data = @(x, U)obj.dataCache(x, varargin{1}, U);
          obj.state = varargin{2};
        end
        if ~isempty(obj.loc)
          if nargin(obj.loc) > 1 % loc(x,t)
            obj.matrix = [];
            obj.idx = obj.fesTrial.mesh.topology.isBoundary(@(x)obj.loc(x, varargin{1}));
          else
            if strcmp(obj.idx, ':')
              obj.idx = obj.fesTrial.mesh.topology.isBoundary(@(x)obj.loc(x));
            end
          end
        end
      end
    end
  end
  methods % assemble
    function assemble(obj)
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
        I = (r.*c==0); if any(I(:)), r(I) = []; c(I) = []; e(I) = []; end
        obj.matrix = obj.matrix + sparse(r(:), c(:), e(:), M, N);
        if k>1
          fprintf(repmat('\b',1,length(s)));
        end
        s = sprintf('progress assembly LHS: %d / %d', k, nBlock); fprintf(s);          
      end
      fprintf('\n');
    end
    function R = integrate(obj, hasCoeff, basisI, basisJ, k)
      [~, weights] = obj.fesTrial.getQuadData(obj.codim);
      [~,~,jac] = obj.fesTrial.evalTrafoInfo([], obj.codim, {k}); % nExnP
      if hasCoeff
        coeff = obj.fesTrial.evalFunction(obj.data, [], obj.codim, obj.state, {k}); % nExnP
        dX = bsxfun(@times, coeff.*abs(jac), weights'); % nExnP
      else
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
      end
      nE = size(basisI, 1); nBI = size(basisI, 2); nBJ = size(basisJ,2); nP = size(basisI,3);
      basisI = reshape(basisI, nE, nBI, nP, []);
      basisJ = reshape(basisJ, nE, nBJ, nP, []);
      R = sum(bsxfun(@times, permute(basisI, [1 2 5 3 4]), ...
                             permute(basisJ, [1 5 2 3 4])), 5); % nExnBIxnBJxnPx(nC*nD)
      R = sum(bsxfun(@times, R, permute(dX, [1 3 4 2])), 4); % nExnBIxnBJ
    end
  end
end
