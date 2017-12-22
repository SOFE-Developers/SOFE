classdef Functional < SOFE
  properties
    codim
    data, dataCache
    fes
    vector
    loc, idx
    pde
  end
  methods % constructor
    function obj = Functional(data, fes, codim, varargin) % [loc]
      if isnumeric(data) && numel(data)<4
        data = @(x)data+zeros(size(x,1),numel(data)); 
      end
      obj.dataCache = data;
      obj.data = data;
      obj.fes = fes;
      obj.fes.register(obj);
      obj.codim = codim;
      obj.idx = ':';
      if nargin > 3
        obj.loc = varargin{1};
      elseif codim == 1
        obj.loc = @(x)~obj.fes.fixB(x);
      end
      obj.pde.state = []; obj.pde.dState = [];
    end
    function notify(obj, varargin) % [time, state, dState]
      if nargin < 2
        obj.vector = [];
        obj.idx = ':';
      else
        if ~isnumeric(obj.dataCache)
          if nargin(obj.dataCache) == 2 % f(x,t)
            obj.vector = [];
            obj.data = @(x)obj.dataCache(x, varargin{1});  
          elseif nargin(obj.dataCache) == 3 % f(x,t,U)
            obj.vector = [];
            obj.data = @(x, U)obj.dataCache(x, varargin{1}, U);
          elseif nargin(obj.dataCache) == 4 % f(x,t,U,d)
            obj.vector = [];
            obj.data = @(x, U, d)obj.dataCache(x, varargin{1}, U, d);
          end
        end
        if nargin > 2
          obj.pde.state = varargin{2};
        end
        if nargin > 3
          obj.pde.dState = varargin{3};
        end
        if ~isempty(obj.loc)
          if nargin(obj.loc) > 1 % loc(x,t)
            obj.vector = [];
            obj.idx = obj.fes.mesh.isBoundary(@(x)obj.loc(x, varargin{1}));
          else
            if strcmp(obj.idx, ':')
              obj.idx = obj.fes.mesh.isBoundary(@(x)obj.loc(x));
            end
          end
        end
      end
    end
  end
  methods
    function assemble(obj)
      % single accumulation
      if ~isempty(obj.vector), return, end
      obj.vector = zeros(obj.fes.getNDoF(), 1);
      if ~any(obj.idx), return, end
      nBlock = obj.fes.nBlock(obj.codim+1);
      re = cell(nBlock,1);
      for k = 1:nBlock
        I = obj.fes.getBlock(obj.codim, k);
        e = []; r = [];
        if ~isempty(I)
          e = obj.assembleOp(k);
          r = abs(obj.fes.getDoFMap(obj.codim, I))';
          if ~ischar(obj.idx)
            e = e(obj.idx(I),:); r = r(obj.idx(I),:);
          end
        end
        I = (r==0); if any(I(:)), r(I) = []; e(I) = []; end %#ok<AGROW>
        re{k} = [r(:), e(:)];
        if k>1
          if k>2
            fprintf(repmat('\b',1,length(s)));
          end
          s = sprintf('progress assembly RHS: %d / %d', k, nBlock);fprintf(s);          
        end
      end
      if k>1, fprintf('\n'); end
      re = cell2mat(re);
      obj.vector = accumarray(re(:,1), re(:,2), size(obj.vector));
    end
    function assemble_(obj)
      % incremental accumulation
      if ~isempty(obj.vector), return, end
      obj.vector = zeros(obj.fes.getNDoF(), 1);
      if ~any(obj.idx), return, end
      nBlock = obj.fes.nBlock(obj.codim+1);
      for k = 1:nBlock
        I = obj.fes.getBlock(obj.codim, k);
        e = []; r = [];
        if ~isempty(I)
          e = obj.assembleOp(k);
          r = abs(obj.fes.getDoFMap(obj.codim, I))';
          if ~ischar(obj.idx)
            e = e(obj.idx(I),:); r = r(obj.idx(I),:);
          end
        end
        I = (r==0); if any(I(:)), r(I) = []; e(I) = []; end %#ok<AGROW>
        obj.vector = obj.vector + accumarray(r(:), e(:), size(obj.vector));
        if k>1
          if k>2
            fprintf(repmat('\b',1,length(s)));
          end
          s = sprintf('progress assembly RHS: %d / %d', k, nBlock);fprintf(s);          
        end
      end
      if k>1, fprintf('\n'); end
    end
  end
end
