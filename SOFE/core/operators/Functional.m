classdef Functional < SOFE
  properties
    codim
    data, dataCache
    fes
    matrix
    loc
    hasCoeff = true;
    matrixFree = false
    preMatrix % row - column - entry
  end
  methods % constructor
    function obj = Functional(data, fes, codim, varargin) % [loc]
      if isnumeric(data) && size(data,1)==1
        obj.dataCache = @(x)data+zeros(size(x,1),numel(data)); 
      else
        obj.dataCache = data;
      end
      obj.data = data;
      obj.fes = fes;
      obj.fes.register(obj);
      obj.codim = codim;
      if ~isempty(varargin)
        obj.loc = varargin{1};
      elseif codim == 1
        obj.loc = @(x)~obj.fes.fixB(x);
      end
    end
    function notify(obj, message, varargin) % [time]
      switch message
        case 'fesChanged'
          if ~obj.matrixFree
            obj.matrix = [];
          else
            % uniform refine (so far, TODO: adaptive)
            nE = obj.fes.mesh.topology.getNumber('0');
            N = nE/size(obj.preMatrix{2},1);
            switch obj.fes.element.conformity
              case {'HDiv','HRot'}
                scal = 1/sqrt(N);
              otherwise
                scal = 1/N;
            end
            obj.preMatrix{2} = repmat(scal*obj.preMatrix{2}, N, 1); % nExnB
            coeff = obj.dataCache(obj.fes.mesh.getCenter('0'));
            e = obj.preMatrix{2}.*repmat(coeff,1,size(obj.preMatrix{2},2)/size(coeff,2)); % nExnB
            r = obj.fes.getDoFMap(0)'; % nExnB
            e(r<0) = -e(r<0);
            obj.matrix = accumarray(abs(r(:)), e(:));
          end
          obj.notifyObservers('fcChanged');
        case 'stateChanged'
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
        otherwise
          error('Unknown message');
      end
    end
  end
  methods % assemble & integrate
    function assemble(obj)
      if ~isempty(obj.matrix), return, end
      obj.matrix = zeros(obj.fes.getNDoF(), 1);
      if ~isempty(obj.loc)
        idx = obj.fes.mesh.isBoundary(@(x)obj.loc(x));
      else
        idx = ':';
      end
      if ~any(idx), return, end
      nBlock = obj.fes.nBlock(obj.codim+1);
      re = cell(nBlock,1);
      for k = 1:nBlock
        I = obj.fes.getBlock(obj.codim, k);
        e = []; r = [];
        if ~isempty(I)
          e = obj.assembleOp(k);
          r = abs(obj.fes.getDoFMap(obj.codim, {k}))';
          if ~strcmp(idx,':')
            e = e(idx(I),:);
            r = r(idx(I),:);
          end
        end
        I = (r==0); if any(I(:)), r(I) = []; e(I) = []; end %#ok<AGROW>
        if obj.matrixFree
          assert(k==1, 'No blocking for matrix free coarse grid');
          obj.preMatrix = {r; e.*sign(obj.fes.getDoFMap(obj.codim, {k}))'};
        end
        re{k} = [r(:), e(:)];
        if k>1
          if k>2
            fprintf(repmat('\b',1,length(s)));
          end
          s = sprintf(['progress assembly ', class(obj), ': %d / %d'], k, nBlock);
          fprintf(s);
        end
      end
      if k>1, fprintf('\n'); end
      re = cell2mat(re);
      obj.matrix = accumarray(re(:,1), re(:,2), size(obj.matrix));
    end
    function R = integrate(obj, basis, k)
      [~, weights] = obj.fes.element.getQuadData(obj.codim); % nPx1
      if obj.hasCoeff
        if isnumeric(obj.data)
          if size(obj.data,1)==1
            R = permute(obj.data, [3 1 2]); % 1x1xnC
          else
            R = obj.fes.evalDoFVector(obj.data, [], obj.codim, 0, {k}); % nExnPxnC
          end
        else
          try S = obj.observers{1}.evalState(k); catch, S = []; end
          R = obj.fes.evalFunction(obj.data, [], obj.codim, S, {k}); % nExnPxnC
        end
      else
        R = 1.0;
      end
      if ~isempty(basis)
        [~,~,jac] = obj.fes.evalTrafoInfo([], obj.codim, {k}); % nExnP
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
        R = sum(bsxfun(@times, permute(R, [1 4 2 3]), basis), 4); % nExnBxnP
        R = sum(bsxfun(@times, R, permute(dX, [1 3 2])), 3); % nExnB
      else
        R = repmat(R,numel(obj.fes.getBlock(k)),1,1);
      end  
    end
  end
end
