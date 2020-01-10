classdef Functional < SOFE
  properties
    codim
    data, dataCache
    fes
    matrix
    loc
    hasCoeff = true;
    matrixFree = false
    A0
  end
  methods % constructor
    function obj = Functional(data, fes, codim, varargin) % [loc]
      obj.setData(data);
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
            if isempty(obj.A0), return; end
            dim = obj.fes.element.dimension;
            hFactor = 0.5;
            switch obj.fes.element.conformity
              case {'HRot'}
                scal = hFactor^(dim-1);
                if obj.fes.element.isSimplex()
                  assert(dim==2, 'TODO');
                  obj.A0 = [repmat(scal*obj.A0, 3, 1, 1); -scal*obj.A0]; % nExnBxnC
                else
                  obj.A0 = repmat(scal*obj.A0, 2^dim, 1, 1); % nExnBxnC
                end
              case {'HDiv'}
                scal = hFactor;
                if obj.fes.element.isSimplex()
                  assert(dim==2, 'TODO');
                  obj.A0 = [repmat(scal*obj.A0, 3, 1, 1); -scal*obj.A0]; % nExnBxnC
                else
                  obj.A0 = repmat(scal*obj.A0, 2^dim, 1, 1); % nExnBxnC
                end
              otherwise
                scal = hFactor^dim;
                obj.A0 = repmat(scal*obj.A0, 2^dim, 1, 1); % nExnBxnC
            end
            %
            coeff = obj.dataCache(obj.fes.mesh.getCenter('0')); % nExnC
            e = sum(obj.A0.*permute(coeff, [1 3 2]), 3); % nExnB
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
    function setData(obj, data)
      if isnumeric(data) && size(data,1)==1
        obj.dataCache = @(x)data+zeros(size(x,1),numel(data)); 
      else
        obj.dataCache = data;
      end
      obj.data = data;
      obj.matrix = [];
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
      obj.A0 = cell(nBlock,1);
      for k = 1:nBlock
        I = obj.fes.getBlock2(obj.codim, k);
        e = []; r = [];
        if ~isempty(I)
          e = obj.assembleOp(k);
          r = abs(obj.fes.getDoFMap(obj.codim, {k}))';
          if ~strcmp(idx,':')
            e = e(idx(I(1):I(2)),:);
            r = r(idx(I(1):I(2)),:);
          end
        end
        if obj.matrixFree
          obj.A0{k} = e.*sign(obj.fes.getDoFMap(obj.codim, {k}))'; % nExnB
        else
          I = (r==0); if any(I(:)), r(I) = []; e(I) = []; end %#ok<AGROW>
          re{k} = [r(:), e(:)];
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
      if obj.matrixFree
        obj.A0 = cell2mat(obj.A0);
      else
        re = cell2mat(re);
        obj.matrix = accumarray(re(:,1), re(:,2), size(obj.matrix));
      end
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
        if obj.matrixFree
          R = bsxfun(@times, permute(R, [1 4 2 3]), basis); % nExnBxnPxnC
        else
          R = sum(bsxfun(@times, permute(R, [1 4 2 3]), basis), 4); % nExnBxnP
        end
        R = sum(permute(bsxfun(@times, R, permute(dX, [1 3 2])), [1 2 4 3]), 4); % nExnB[xnC]
      else
        nE = obj.fes.getBlock2(k);
        R = repmat(R,nE(2)-nE(1)+1,1,1);
      end  
    end
  end
end
