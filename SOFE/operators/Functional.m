classdef Functional < SOFEClass
  properties
    codim
    data, dataCache
    feSpace
    vector
    state
    loc, idx
  end
  methods % constructor
    function obj = Functional(data, feSpace, codim, varargin) % [loc]
      if isreal(data), data = @(x)data+0*x(:,1); end
      obj.dataCache = data;
      obj.data = data;
      obj.feSpace = feSpace;
      obj.feSpace.register(obj);
      obj.codim = codim;
      obj.idx = ':';
      if nargin > 3
        obj.loc = varargin{1};
      elseif codim == 1
        obj.loc = @(x)~obj.feSpace.fixB(x);
      end
    end
    function notify(obj, varargin) % [time, state]
      if nargin < 2
        obj.vector = [];
        obj.idx = ':';
      else
        if nargin(obj.dataCache) == 2 % f(x,t)
          obj.vector = [];
          obj.data = @(x)obj.dataCache(x, varargin{1});  
        elseif nargin(obj.dataCache) == 3 % f(x,t,U)
          obj.vector = [];
          obj.data = @(x, U)obj.dataCache(x, varargin{1}, U);
          obj.state = varargin{2};
        end
        if ~isempty(obj.loc)
          if nargin(obj.loc) > 1 % time dependent location
            obj.vector = [];
            obj.idx = obj.feSpace.mesh.topology.isBoundary(@(x)obj.loc(x, varargin{1}));
          else
            if strcmp(obj.idx, ':')
              obj.idx = obj.feSpace.mesh.topology.isBoundary(@(x)obj.loc(x));
            end
          end
        end
      end
    end
  end
  methods
    function assemble(obj)
      if ~isempty(obj.vector), return, end
      M = obj.feSpace.getNDoF();
      obj.vector = zeros(M, 1);
      if ~any(obj.idx), return, end
      nBlock = obj.feSpace.mesh.nBlock;
      for k = 1:nBlock
        [r,e] = obj.assembleBlock(k);
        if ~isempty(r) && numel(obj.idx) > 1
          I = obj.feSpace.mesh.getBlock(obj.codim, k);
          r = r(obj.idx(I(1):I(2)),:);
          e = e(obj.idx(I(1):I(2)),:);
        end
        I = r==0;
        if any(I(:))
          r(I) = []; e(I) = [];
        end
        obj.vector = obj.vector + accumarray(r(:), e(:), [M, 1]);
        if nBlock > 1
          if k>1
            fprintf(repmat('\b',1,length(s)));
          end
          s=sprintf('progress assembly RHS: %d / %d', k, nBlock);fprintf(s);          
        end
      end
      if nBlock > 1, fprintf('\n'); end
    end
    function [r,e] = assembleBlock(obj, k)
      I = obj.feSpace.mesh.getBlock(obj.codim, k);
      if ~isempty(I)
        e = obj.assembleOp(k);
        r = abs(obj.feSpace.getDoFMap(obj.codim, (I(1):I(2))'))';
      else
        r = []; e = [];
      end
    end
  end
end
