classdef MeshTopology < SOFEClass
  properties
    dimW
    dimP
    nodes
    connectivity
    globalSearcher
    %
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
      try
        obj.globalSearcher = GlobalSearcher(obj, 10);
      catch
        warning('Building GlobalSearcher failed');
      end
    end
  end
  methods % obj is observed
    function register(obj, observer)
      obj.observers = [obj.observers, {observer}];
    end
    function notifyObservers(obj)
      try
        obj.globalSearcher = GlobalSearcher(obj, 10);
      catch
        warning('Building GlobalSearcher failed');
      end
      for i = 1:numel(obj.observers)
        obj.observers{i}.notify();
      end
    end
  end
  methods % mesh information
    function R = getEntity(obj, dim, varargin)
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
        R = reshape(obj.nodes(obj.getEntity(dim),:), obj.getNumber(dim), [], size(obj.nodes,2)); % nExnVxnW
        R = loc(permute(mean(R,2), [1 3 2])); % nE
      else
        R = true(obj.getNumber(dim), 1);
      end
    end
    function R = getNumber(obj, dim)
      R = size(obj.getEntity(dim), 1);
    end
    function R = getNode2Elem(obj)
      nE = obj.getNumber(obj.dimP);
      elem = obj.getEntity(obj.dimP);
      elemIndex = repmat((1:nE)', 1, size(elem,2));
      elem = elem(:);
      [dmy,I] = sort(elem);
      count = accumarray(elem,1);
      maxCount = max(count);
      upperTri = triu(repmat((1:maxCount)',1,maxCount));
      count = upperTri(:,count); count = count(:);
      count(count==0) = [];
      R = accumarray([elem(I), count], elemIndex(I));
    end
    function R = getCenter(obj, dim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      if dim == 0
        vertices = obj.getEntity(0);
        R = vertices(I,:);
        return;
      end
      R = obj.getEntity(dim); R = R(I,:);
      [nE, nV] = size(R);
      vertices = obj.getEntity(0);
      R = permute(sum(reshape(vertices(R,:), nE, nV, []),2)/nV,[1 3 2]);
    end
    function R = isBoundary(obj, varargin)
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
    function R = isSurface(obj, varargin)
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
  methods(Static = true)
    function R = getTopology(nodes, elem, dimP)
      switch size(elem, 2)
        case 2
          R = MeshTopologyInt(nodes, elem, dimP);
        case 3
          R = MeshTopologyTri(nodes, elem, dimP);
        case 4
          if dimP == 2
            R = MeshTopologyQuad(nodes, elem, dimP);
          else
            R = MeshTopologyTet(nodes, elem, dimP);
          end
        case 8
          R = MeshTopologyHex(nodes, elem, dimP);
      end
    end
    function R = getShapeElement(N, dimP)
      switch N
        case 2
          R = P1(1);
        case 3
          R = P1(2);
        case 4
          if dimP == 2
            R = Q1(2);
          else
            R = P1(3);
          end
        case 8
          R = Q1(3);
      end
    end
  end
end