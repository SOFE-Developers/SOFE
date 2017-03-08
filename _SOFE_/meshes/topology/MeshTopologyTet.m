classdef MeshTopologyTet < MeshTopology
  methods % constructor
    function obj = MeshTopologyTet(nodes, elem, dimP)
      elem = MeshTopologyTet.renumber(nodes, elem);
      obj = obj@MeshTopology(nodes, elem, dimP);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{3,1} = [obj.connectivity{4,1}(:,[1,2,3]); ...
                               obj.connectivity{4,1}(:,[1,2,4]); ...
                               obj.connectivity{4,1}(:,[2,3,4]); ...
                               obj.connectivity{4,1}(:,[1 3 4])];
      [obj.connectivity{3,1}, dmy, e2F] = unique(sort(obj.connectivity{3,1},2),'rows');    
      obj.connectivity{2,1} = [obj.connectivity{4,1}(:,[1,2]); obj.connectivity{4,1}(:,[2,3]); ...
                               obj.connectivity{4,1}(:,[1,3]); obj.connectivity{4,1}(:,[1 4]); ...
                               obj.connectivity{4,1}(:,[2,4]); obj.connectivity{4,1}(:,[3 4])];
      [obj.connectivity{2,1}, dmy, e2Ed] = unique(sort(obj.connectivity{2,1},2),'rows');    
      obj.connectivity{4,3} = reshape(e2F,[], 4);
      obj.connectivity{4,2} = reshape(e2Ed,[], 6);
      obj.connectivity{3,2} = obj.getFace2Edge();
      %
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
      obj.connectivity{4,4} = (1:size(obj.connectivity{4,1},1))';
    end
  end
  methods % connectivity information
    function R = getElem2Face(obj)
      R = obj.connectivity{4,3};
    end
    function R = getElem2Edge(obj)
      R = obj.connectivity{4,2};
    end
    function R = getFace2Edge(obj)
      R = zeros(obj.getNumber(2), 3);
      e2F = obj.getElem2Face();
      e2E = obj.getElem2Edge();
      orientF = obj.getOrientation(3,2);
      orient1 = orientF(:,:,1);
      orient2 = orientF(:,:,2);
      [dmy, ind] = unique(e2F);
      orient1 = orient1(ind); orient2 = orient2(ind);
      [elems, type] = ind2sub([obj.getNumber(3), 4], ind);
      nodeIxAtFace = [1 2 3; 1 5 4; 2 6 5; 3 6 4];
      for t = 1:4
        for k = 0:2
          I = (type == t) & (orient1 == k+1);
          R(I,:) = e2E(elems(I), circshift(nodeIxAtFace(t,:)',-k));
        end
      end
      I = (orient2 == -1);
      R(I,:) = R(I, [3 2 1]);
    end
    function R = getFace2Elem(obj)
      nE = obj.getNumber(3); nF = obj.getNumber(2);
      orient = obj.getNormalOrientation();
      R = full(sparse(obj.getElem2Face(), 0.5*(3-orient), repmat((1:nE)',1,4),nF,2));
    end
    function R = getOrientation(obj, dim1, dim2)
      R = [];
      switch dim2
        case 1
          switch dim1
            case 3
              e = obj.getEntity(3);
              R = ones(size(e,1), 6);
              R(e(:,1)>e(:,2),1) = -1;
              R(e(:,2)>e(:,3),2) = -1;
              R(e(:,1)>e(:,3),3) = -1;
              R(e(:,1)>e(:,4),4) = -1;
              R(e(:,2)>e(:,4),5) = -1;
              R(e(:,3)>e(:,4),6) = -1;
            case 2
              R = ones(obj.getNumber(2), 3);
          end
        case 2
          e = obj.getEntity(3);
          R = ones(2, size(e,1), 4); % ! two orient flags
          face = [1 2 3; 1 2 4; 2 3 4; 1 3 4];
          for i = 1:4
              [~, R(1,:,i)] = min(e(:,face(i,:)),[],2);
              [~, P] = sort(e(:,face(i,:)),2);
              even = (P(:,1)<P(:,2) & P(:,2)<P(:,3)) | ...
                     (P(:,2)<P(:,3) & P(:,3)<P(:,1)) | ...
                     (P(:,3)<P(:,1) & P(:,1)<P(:,2));
              R(2,~even,i) = -1;
          end
          R = permute(R,[2 3 1]); % nExnFx2
      end
    end
    function R = getNormalOrientation(obj)
      orient = obj.getOrientation([], 2); % nEx4x2
      R = ones(obj.getNumber(3), 4); % nExnF
      R(:,[2 3]) =  orient(:,[2 3],2);
      R(:,[1 4]) = -orient(:,[1 4],2);
    end
  end
  methods % mesh information
    function R = getMeasure(obj, dim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      ee = obj.getEntity(dim); ee = ee(I,:);
      nodes = obj.getEntity(0);
      switch dim
        case 3 % element
          v1 = nodes(ee(:,2),:) - nodes(ee(:,1),:);
          v2 = nodes(ee(:,3),:) - nodes(ee(:,1),:);
          v3 = nodes(ee(:,4),:) - nodes(ee(:,1),:);
          R = ((v1(:,1).*v2(:,2).*v3(:,3) + v1(:,2).*v2(:,3).*v3(:,1)+v1(:,3).*v2(:,1).*v3(:,2)) ...
          - (v1(:,3).*v2(:,2).*v3(:,1)+v1(:,2).*v2(:,1).*v3(:,3)+v1(:,1).*v2(:,3).*v3(:,2)))/6;
        case 2 % face
          v1 = nodes(ee(:,2),:) - nodes(ee(:,1),:);
          v2 = nodes(ee(:,3),:) - nodes(ee(:,1),:);
          R = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))/2;
        case 1 % edge
          v = nodes(ee(:,2),:) - nodes(ee(:,1),:);
          R = sum(v.^2,2).^0.5;
       end
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(3); nN = obj.getNumber(0);
      obj.nodes = [obj.nodes; obj.getCenter(1)];
      newNodeNumber = (nN+1 : nN+obj.getNumber(1));
      el = [el newNodeNumber(obj.connectivity{4,2})];
      obj.connectivity{4,1} = [el(:,[1 5 7 8]); el(:,[2 6 5 9]); ...
                               el(:,[3 7 6 10]); el(:,[4 9 8 10]); ...
                               el(:,[5 6 7 9]); el(:,[5 7 8 9]); ...
                               el(:,[7 8 9 10]); el(:,[7 9 6 10])];
      obj.updateConnectivity();
      obj.notifyObservers();
    end
  end
  methods % display
    function show(obj, varargin)
      c = caxis();
      fc = obj.getEntity(2);
      I = obj.isSurface(varargin{:});
      h = trimesh(fc(I,:), obj.nodes(:,1), obj.nodes(:,2), obj.nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
      axis equal, axis tight, caxis(c);
    end
    function showInner(obj)
      nodes = obj.getEntity(0);
      h = trimesh(obj.getEntity(2),nodes(:,1),nodes(:,2),nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
    end
    function showTets(obj, I, flag)
      e2F = obj.getElem2Face();
      e2Ed = obj.getElem2Edge();
      elem = obj.getEntity(3);
      II = cell(4,1);
      II{4} = I;
      II{3} = unique(e2F(I,:));
      II{2} = unique(e2Ed(I,:));
      II{1} = unique(elem(I,:));
      face = obj.getEntity(2);
      nodes = obj.getEntity(0);
      h = trimesh(face(II{3},:),nodes(:,1),nodes(:,2),nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
      if flag
        for i = 1:numel(flag)
          dim = str2double(flag(i));
          obj.showEntity(dim, II{dim+1});
        end
      end
    end
    function tetramesh(obj, varargin)
      if nargin>1, I = varargin{1}; else I=':'; end
      tet = obj.getEntity(3);
      tetramesh(tet(I,:), obj.getEntity(0));
    end
    function showEntity(obj, dim, varargin) % for debug
      if nargin < 3, I = (1:obj.getNumber(dim))'; else I = varargin{1}; end
      center = obj.getCenter(dim);
      switch dim
        case 3
          color = [1 0 1];
        case 2
          color = [0 0 0];
        case 1
          color = [0 0.5 1];
        case 0
          color = [0 1 0.5];
      end
      text(center(I,1), center(I,2), center(I,3), num2str(I(:)),'Color',color,'FontSize', 18);
    end
    function showNodeVector(obj, U, varargin)
      fc = obj.getEntity(2);
      if nargin > 2
        I = obj.isSurface(varargin{1});
      else
        I = obj.isBoundary();
      end
      nodes = obj.getEntity(0);
      h = trimesh(fc(I,:),nodes(:,1),nodes(:,2),nodes(:,3), U(1:obj.getNumber(0)));
      set(h,'facecolor','interp','edgecolor','k');
      axis equal, axis tight
    end
  end
  methods(Static = true)
    function R = getQuadRule(order)
      R{4} = GaussPoint();
      R{3} = GaussInt(order);
      R{2} = GaussTri(order);
      R{1} = GaussTet(order);
    end
    function R = getBarycenterRef()
      R = [1 1 1]/4;
    end
    function R = isFeasible(points)
      tol = 1e-12;
      R = (all(points>-tol, 2) & 1-sum(points,2)>-tol);
    end
    function R = getEquiPoints(order)
      ls = linspace(0,1, order+1);
      [X, Y, Z] = meshgrid(ls, ls, ls);
      R = [X(:) Y(:) Z(:)];
      R = R(sum(R,2)<=1,:);
    end
    function R = renumber(nodes, elem)
      v1 = nodes(elem(:,2),:) - nodes(elem(:,1),:);
      v2 = nodes(elem(:,3),:) - nodes(elem(:,1),:);
      v3 = nodes(elem(:,4),:) - nodes(elem(:,1),:);
      I = ((v1(:,1).*v2(:,2).*v3(:,3) + v1(:,2).*v2(:,3).*v3(:,1)+v1(:,3).*v2(:,1).*v3(:,2)) ...
      - (v1(:,3).*v2(:,2).*v3(:,1)+v1(:,2).*v2(:,1).*v3(:,3)+v1(:,1).*v2(:,3).*v3(:,2)))/6;
      I = I<0;
      if any(I)
        fprintf('Elements renumbered!\n');
        elem(I,:) = elem(I, [2 1 3 4]);
      end
      R = elem;
    end
  end
end