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
    function obj = MeshTopology(nodes, elem, dimP)
      obj.nodes = nodes;
      obj.dimW = size(nodes,2);
      obj.dimP = dimP;
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.observers = {};
      obj.buildGlobalSearcher();
    end
    function buildGlobalSearcher(obj)
      try
        obj.globalSearcher = GlobalSearcher(obj);
      catch
        fprintf('Building GlobalSearcher failed\n');
      end
    end
  end
  methods % obj is observed
    function register(obj, observer)
      obj.observers = [obj.observers, {observer}];
    end
    function notifyObservers(obj)
      obj.buildGlobalSearcher();
      for i = 1:numel(obj.observers)
        obj.observers{i}.notify();
      end
    end
  end
  methods % mesh information
    function R = getEntity(obj, dim, varargin) % [I]
      I = ':'; if nargin > 2, I = varargin{1}; end
      if dim == 0
        R = obj.nodes(I,:);
      else
        R = obj.connectivity{dim+1}(I,:);
      end
    end
    function R = findEntity(obj, dim, varargin) % [loc]
      if nargin > 2
        loc = varargin{1};
        R = any(reshape(loc(obj.nodes(obj.getEntity(dim),:)), obj.getNumber(dim), []), 2);
      else
        R = true(obj.getNumber(dim), 1);
      end
    end
    function R = findEntityC(obj, dim, varargin) % [loc]
      if nargin > 2
        loc = varargin{1};        
        R = loc(obj.getCenter(dim)); % nE
      else
        R = true(obj.getNumber(dim), 1);
      end
    end
    function R = getNumber(obj, dim, varargin) % [I]
      R = size(obj.getEntity(dim), 1, varargin{:});
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
    function R = getMeasure(obj, dim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      ee = obj.getEntity(dim); ee = ee(I,:);
      switch dim
        case 3
          v1 = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          v2 = obj.nodes(ee(:,3),:) - obj.nodes(ee(:,1),:);
          v3 = obj.nodes(ee(:,4+(size(ee,2)==8)),:) - obj.nodes(ee(:,1),:);
          R = ((v1(:,1).*v2(:,2).*v3(:,3) + v1(:,2).*v2(:,3).*v3(:,1)+v1(:,3).*v2(:,1).*v3(:,2)) ...
          - (v1(:,3).*v2(:,2).*v3(:,1)+v1(:,2).*v2(:,1).*v3(:,3)+v1(:,1).*v2(:,3).*v3(:,2)));
          if size(ee,2)==4
            R = R/6;  
          else
            if all(all(obj.nodes(ee(:,1),:) + v1+v2+v3 - obj.nodes(ee(:,8),:)>1e-12))
              warning('! Volume only valid for parallelepiped !');
            end
          end
        case 2
          v1 = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          v2 = obj.nodes(ee(:,3),:) - obj.nodes(ee(:,1),:);
          R = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1));
          if size(ee,2)==3
            R = R/2;
          else
            if all(all(obj.nodes(ee(:,1),:) + v1+v2 - obj.nodes(ee(:,4),:)>1e-12))
              warning('! Area only valid for parallelograms !');
            end
          end
        case 1
          v = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          R = sum(v.^2,2).^0.5;
       end
    end
    function R = getDiam(obj)
      R = [min(obj.nodes); max(obj.nodes)];
    end
    function R = isBoundary(obj, varargin) % [loc]
      e2F = obj.getElem2Face(); % nExnF
      R = accumarray(e2F(e2F>0),1, [obj.getNumber(obj.dimP-1) 1])==1; % nFx1
      if nargin > 1
        if ~isempty(varargin{1})
          R = find(R); % nBFx1
          I = varargin{1}(obj.getCenter(obj.dimP-1, R)); nC = size(I,2); % nBFxnC
          R = repmat(R, 1, nC); % nBFxnC
          col = ones(size(R,1),1)*[1 2 3]; % nBFxnC
          R = accumarray([R(I), col(I)], 1, [obj.getNumber(obj.dimP-1), nC])>0; % nFxnC
        else
          R = [];
          return
        end
      end
    end
    function R = isSurface(obj, varargin) % [loc]
      if nargin > 1 && ~ischar(varargin{1})
        idx = find(varargin{:}(obj.getEntity(0)));
      else
        idx = 1:obj.getNumber(0);
      end
      goodElem = any(ismember(obj.getEntity(obj.dimP), idx),2); % feasible elements
      E2F = obj.getElem2Face();
      E2F = E2F(goodElem,:);
      uE2F = unique(E2F(:));
      I = hist(E2F(:),uE2F)==1;
      R = full(sparse(uE2F(I), 1, true, obj.getNumber(obj.dimP-1), 1));
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
      nE = obj.getNumber(obj.dimP); nF = obj.getNumber(obj.dimP-1);
      orient = 0.5*(3-obj.getNormalOrientation());
      R = full(sparse(obj.getElem2Face(), orient, repmat((1:nE)',1,size(orient,2)),nF,2));
      type = full(sparse(obj.getElem2Face(), orient, ones(nE,1)*(1:size(orient,2)),nF,2));
    end
    function R = getNodePatch(obj, dim)
      if dim>0
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
        R = accumarray([entity(I), count], idx(I));
      else
        segm = obj.getEntity(1);
        if dim < 0 % on boundary
          iB = obj.isBoundary;
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
  methods % mesh operations
    function rotate(obj, alpha)
      if obj.dimW ~= 2, error('!Dimension must be 2!'); end
      A = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
      obj.applyLinearMap(A);
    end
    function scale(obj, a)
      A = diag(a)*eye(obj.dimW);
      obj.applyLinearMap(A);
    end
    function shift(obj, vec)
      obj.nodes = bsxfun(@plus, obj.nodes, vec(:)');
      obj.notifyObservers();
    end
    function applyLinearMap(obj, A)
      obj.nodes = (A*obj.nodes')';
      obj.notifyObservers();
    end
  end
end