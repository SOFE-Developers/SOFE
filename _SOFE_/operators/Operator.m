classdef Operator < SOFEClass
  properties
    codim = 0;
    data, dataCache
    feSpaceTrial, feSpaceTest
    isGalerkin
    matrix
    state
    loc, idx
  end
  methods % constructor
    function obj = Operator(data, feSpaceTrial, varargin) % [feSpaceTest loc]
      if isreal(data), data = @(x)data+0*x(:,1); end
      obj.dataCache = data;
      obj.data = data;
      obj.feSpaceTrial = feSpaceTrial;
      obj.feSpaceTrial.register(obj);
      if nargin > 2 && ~isempty(varargin{1})
        obj.feSpaceTest = varargin{1};
        obj.feSpaceTest.register(obj);
        obj.isGalerkin = false;
        % sync quadRules
        if obj.feSpaceTrial.quadRule{1}.order > obj.feSpaceTest.quadRule{1}.order
          obj.feSpaceTest.setQuadRule(obj.feSpaceTrial.quadRule);
        else
          obj.feSpaceTrial.setQuadRule(obj.feSpaceTest.quadRule);
        end
      else
        obj.feSpaceTest = feSpaceTrial;
        obj.isGalerkin = true;
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
          obj.data = @(x)obj.dataCache(x, varargin{1});
          obj.matrix = [];
        elseif nargin(obj.dataCache) == 3 % f(x,t,U)
          obj.data = @(x, U)obj.dataCache(x, varargin{1}, U);
          obj.state = varargin{2};
          obj.matrix = [];
        end
        if ~isempty(obj.loc)
          if nargin(obj.loc) > 1 % loc(x,t)
            obj.matrix = [];
            obj.idx = obj.feSpaceTrial.mesh.topology.isBoundary(@(x)obj.loc(x, varargin{1}));
          else
            if strcmp(obj.idx, ':')
              obj.idx = obj.feSpaceTrial.mesh.topology.isBoundary(@(x)obj.loc(x));
            end
          end
        end
      end
    end
  end
  methods % assemble
    function assemble(obj)
      if ~isempty(obj.matrix), return, end
      M = obj.feSpaceTest.getNDoF();
      N = obj.feSpaceTrial.getNDoF();
      obj.matrix = sparse(M, N);
      if ~any(obj.idx), return, end
      nBlock = obj.feSpaceTrial.mesh.nBlock;
      for k = 1:nBlock
        [r,c,e] = obj.assembleBlock(k);
        if ~isempty(r) && ~ischar(obj.idx) % numel(obj.idx) > 1
          I = obj.feSpaceTrial.mesh.getBlock(obj.codim, k);
          r = r(obj.idx(I(1):I(2)),:,:);
          c = c(obj.idx(I(1):I(2)),:,:);
          e = e(obj.idx(I(1):I(2)),:,:);
        end
        I = r.*c==0;
        if any(I(:))
          r(I) = []; c(I) = []; e(I) = [];
        end
        obj.matrix = obj.matrix + sparse(r(:), c(:), e(:), M, N);
        if nBlock > 1
          if k>1
            fprintf(repmat('\b',1,length(s)));
          end
          s = sprintf('progress assembly LHS: %d / %d', k, nBlock); fprintf(s);          
        end
      end
      if nBlock > 1, fprintf('\n'); end
    end
    function [r,c,e] = assembleBlock(obj, k)
      I = obj.feSpaceTrial.getBlock(obj.codim, k);
      if ~isempty(I)
        e = obj.assembleOp(k); % nExnBxnB
        c = obj.feSpaceTrial.getDoFMap(obj.codim, (I(1):I(2))'); % nBxnE
        r = obj.feSpaceTest.getDoFMap(obj.codim, (I(1):I(2))'); % nBxnE
        r = repmat(abs(r)',[1 1 size(c,1)]); % nExnBxnB
        c = permute(repmat(abs(c)',[1 1 size(r,2)]), [1 3 2]); % nExnBxnB
      else
        r = []; c = []; e = [];
      end
    end
    function R = integrate(obj, hasCoeff, basisI, basisJ, k)
      [~, weights] = obj.feSpaceTrial.getQuadData(obj.codim);
      [~,~,jac] = obj.feSpaceTrial.evalTrafoInfo([], obj.codim, k); % nExnP
      if hasCoeff
        coeff = obj.feSpaceTrial.evalFunction(obj.data, [], obj.codim, obj.state, k); % nExnP
        dX = bsxfun(@times, coeff.*abs(jac), weights'); % nExnP
      else
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
      end
      nE = size(basisI, 1); nBI = size(basisI, 2); nBJ = size(basisJ,2); nP = size(basisI,3);
      basisI = reshape(basisI, nE, nBI, nP, []);
      basisJ = reshape(basisJ, nE, nBJ, nP, []);
      R = sum(bsxfun(@times, permute(basisI, [1 2 5 3 4]), ...
                             permute(basisJ, [1 5 2 3 4])), 5); % nExnBIxnBJxnP
      R = sum(bsxfun(@times, R, permute(dX, [1 3 4 2])), 4); % nExnBIxnBJ
    end
  end
end
