classdef MeshTopology < SOFE
  properties
    dimW
    dimP
    nodes
    connectivity
    globalSearcher
    observers
  end
  methods % constructor
    function obj = MeshTopology(nodes, dimP)
      obj.nodes = nodes;
      obj.dimW = size(nodes,2);
      obj.dimP = dimP;
      obj.observers = {};
    end
    function R = getGlobalSearcher(obj)
      if isempty(obj.globalSearcher)
        obj.globalSearcher = GlobalSearcher(obj);
      end
      R = obj.globalSearcher;
    end
  end
  methods % obj is observed
    function register(obj, observer)
      obj.observers = [obj.observers, {observer}];
    end
    function notifyObservers(obj)
      obj.globalSearcher = [];
      for i = 1:numel(obj.observers)
        obj.observers{i}.notify();
      end
    end
  end
  methods % mesh information
    function R = getEntity(obj, dim, varargin) % [I]
      I = ':'; if nargin > 2, I = varargin{1}; end
      if ischar(dim), dim  = obj.dimP - str2double(dim); end % dim to codim
      if dim == 0
        R = obj.nodes(I,:);
      else
        R = obj.connectivity{dim+1}(I,:);
      end
    end
    function R = getNumber(obj, dim)
      if ischar(dim), dim  = obj.dimP - str2double(dim); end % dim to codim
      R = numel(obj.connectivity{dim+1, dim+1});
    end
    function R = getCenter(obj, dim, varargin) % [I]
      I = ':'; if nargin > 2, I = varargin{1}; end
      if dim == 0
        R = obj.nodes(I,:);
        return;
      end
      R = obj.getEntity(dim); R = R(I,:);
      [nE, nV] = size(R);
      R = permute(mean(reshape(obj.nodes(R,:), nE, nV, []),2),[1 3 2]);
    end
    function R = isBoundary(obj, varargin) % [loc]
      e2F = obj.getElem2Face(); % nExnF
      R = accumarray(e2F(e2F>0),1, [obj.getNumber('1') 1])==1; % nFx1
      if nargin > 1
        if ~isempty(varargin{1})
          I = varargin{1}(obj.getCenter(obj.dimP-1, R));
          R = repmat(R, 1, size(I,2));
          R(R(:,1)>0,:) = I;
        else
          R = [];
          return
        end
      end
    end
    function R = getBoundary(obj, varargin) % [loc]
      R = obj.isBoundary(varargin{:});
      R = obj.getEntity(1,R);
    end
    function R = isBoundaryNode(obj, varargin) % [loc]
      R = unique(obj.getEntity(1,obj.isBoundary(varargin{:})));
      R = accumarray(R,1,[obj.getNumber(0) 1])>0;
    end
    function R = getBoundaryNode(obj, varargin) % [loc]
      R = obj.isBoundaryNode(varargin{:});
      R = obj.getEntity(0,R);
    end
    function R = isSurface(obj, varargin) % [loc]
      if nargin < 2, R = obj.isBoundary(); return; end
      goodElem = varargin{1}(obj.getCenter(obj.dimP));
      E2F = obj.getElem2Face();
      E2F = E2F(goodElem,:);
      uE2F = unique(E2F(:));
      I = hist(E2F(:),uE2F)==1;
      R = full(sparse(uE2F(I), 1, true, obj.getNumber('1'), 1));
    end
    function R = getSurface(obj, varargin) % [loc]
      R = obj.isSurface(varargin{:});
      R = obj.getEntity(1,R);
    end
  end
  methods % connectivity information
    function R = getElem2Face(obj, varargin)
      R = obj.connectivity{obj.dimP+1,obj.dimP};
      if nargin > 2
        R = R(varargin{1},:);
      end
    end
    function [R, type] = getFace2Elem(obj)
      if isempty(obj.connectivity{obj.dimP,obj.dimP+1})        
        nE = obj.getNumber('0'); nF = obj.getNumber('1');
        orient = 0.5*(3-obj.getNormalOrientation());
        obj.connectivity{obj.dimP,obj.dimP+1}{1} = full(sparse(obj.getElem2Face(), orient, repmat((1:nE)',1,size(orient,2)),nF,2));
        obj.connectivity{obj.dimP,obj.dimP+1}{2} = full(sparse(obj.getElem2Face(), orient, ones(nE,1)*(1:size(orient,2)),nF,2));
      end
      R = obj.connectivity{obj.dimP,obj.dimP+1}{1};
      type = obj.connectivity{obj.dimP,obj.dimP+1}{2};
    end
    function R = getNodePatch(obj, dim)
      if dim>0
        if isempty(obj.connectivity{1,dim+1})
          nE = obj.getNumber(dim);
          entity = obj.getEntity(dim);
          idx = repmat((1:nE)', 1, size(entity,2));
          entity = entity(:);
          [~,I] = sort(entity);
          count = accumarray(entity,1);
          maxCount = max(count);
          upperTri = triu(repmat((1:maxCount)',1,maxCount));
          count = upperTri(:,count); count = count(:);
          count(count==0) = [];
          obj.connectivity{1,dim+1} = accumarray([entity(I), count], idx(I));
        end
        R = obj.connectivity{1,dim+1};
      else
        segm = obj.getEntity(1);
        if dim < 0 % on boundary
          iB = obj.isBoundary();
          segm = segm(iB,:);  
        end
        II = sparse(segm(:,[1 2]),segm(:,[2 1]),segm(:,[2 1]));
        II = leftShiftNonZero(II);
        if dim < 0 % on boundary
          R = full(II(:,1:2));
        else
          R = full(II(:,1:max(sum(II>0,2))));
        end
        R(R(:,1)==0,:) = [];
      end
    end
  end
end