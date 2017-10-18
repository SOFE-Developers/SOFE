classdef MeshTopologyTet < MeshTopology
  methods % constructor
    function obj = MeshTopologyTet(nodes, elem, dimP)
      elem = MeshTopologyTet.renumber(nodes, elem);
      obj = obj@MeshTopology(nodes, dimP);
      obj.updateConnectivity(elem);
    end
    function updateConnectivity(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{3,1} = [elem(:,[1,2,3]); elem(:,[1,2,4]); ...
                               elem(:,[2,3,4]); elem(:,[1 3 4])];
      [obj.connectivity{3,1}, ~, e2F] = unique(sort(obj.connectivity{3,1},2),'rows');    
      obj.connectivity{2,1} = [elem(:,[1,2]); elem(:,[2,3]); elem(:,[1,3]); ...
                               elem(:,[1 4]); elem(:,[2,4]); elem(:,[3 4])];
      [obj.connectivity{2,1}, ~, e2Ed] = unique(sort(obj.connectivity{2,1},2),'rows');    
      obj.connectivity{4,3} = reshape(e2F,[], 4);
      obj.connectivity{4,2} = reshape(e2Ed,[], 6);
      %
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
      obj.connectivity{4,4} = (1:size(obj.connectivity{4,1},1))';
      %
      obj.connectivity{3,2} = obj.getFace2Edge();
    end
  end
  methods % connectivity information
    function R = getElem2Edge(obj)
      R = obj.connectivity{4,2};
    end
    function R = getFace2Edge(obj)
      R = zeros(obj.getNumber(2), 3);
      e2F = obj.getElem2Face();
      e2E = obj.getElem2Edge();
      orientF = obj.getOrientation(3,2);
      [~, ind] = unique(e2F);
      [elems, type] = ind2sub([obj.getNumber(3), 4], ind);
      nodeIxAtFace = [1 2 3; 1 5 4; 2 6 5; 3 6 4];
      for t = 1:4
        for k = 0:2
          I = (type == t) & (abs(orientF(ind)) == k+1);
          R(I,:) = e2E(elems(I), circshift(nodeIxAtFace(t,:)',-k));
        end
      end
      I = orientF(ind)<0;
      R(I,:) = R(I, [3 2 1]);
    end
    function R = getOrientation(obj, dim1, dim2, varargin)
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
            otherwise
              return
          end
        case 2
          switch dim1
            case 3
              e = obj.getEntity(3);
              face = [1 2 3; 1 2 4; 2 3 4; 1 3 4];
              R = zeros(size(e,1),4);
              for i = 1:4
                  [~, R(:,i)] = min(e(:,face(i,:)),[],2);
                  [~, P] = sort(e(:,face(i,:)),2);
                  even = (P(:,1)<P(:,2) & P(:,2)<P(:,3)) | ...
                         (P(:,2)<P(:,3) & P(:,3)<P(:,1)) | ...
                         (P(:,3)<P(:,1) & P(:,1)<P(:,2));
                  R(~even, i) = -R(~even, i);
              end
            otherwise
              return
          end
        otherwise
          return
      end
      if nargin > 3
        R = R(varargin{1},:,:);
      end
    end
    function R = getNormalOrientation(obj)
      R = sign(obj.getOrientation(3, 2)); % nEx4
      R(:,[1 4]) = -R(:,[1 4]);
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(3); nN = obj.getNumber(0);
      obj.nodes = [obj.nodes; obj.getCenter(1)];
      newIndices = (nN+1 : nN+obj.getNumber(1));
      el = [el newIndices(obj.connectivity{4,2})];
      el = [el(:,[1 5 7 8]); el(:,[2 6 5 9]); ...
            el(:,[3 7 6 10]); el(:,[4 9 8 10]); ...
            el(:,[5 6 7 9]); el(:,[5 7 8 9]); ...
            el(:,[7 8 9 10]); el(:,[7 9 6 10])];
      obj.updateConnectivity(el);
    end
  end
  methods % display
    function show(obj, varargin)
      c = caxis();
      fc = obj.getEntity(2);
      Ib = obj.isBoundary(varargin{:});
      Is = obj.isSurface(varargin{:}) & ~Ib;
      h = trimesh(fc(Is,:), obj.nodes(:,1), obj.nodes(:,2), obj.nodes(:,3));
      set(h,'facecolor',[0.5 0.7 0.2],'edgecolor','k');
      hold on
      h = trimesh(fc(Ib,:), obj.nodes(:,1), obj.nodes(:,2), obj.nodes(:,3));
      hold off
      set(h,'facecolor',[0.5 0.8 0.5],'edgecolor','k');
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
      if nargin>1, I = varargin{1}; else, I=':'; end
      tet = obj.getEntity(3);
      tetramesh(tet(I,:), obj.getEntity(0));
    end
    function showEntity(obj, dim, varargin) % for debug
      if nargin < 3, I = (1:obj.getNumber(dim))'; else, I = varargin{1}; end
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
  end
  methods(Static = true)
    function R = getQuadRule(order)
      R{4} = GaussPoint();
      R{3} = GaussInt(order);
      R{2} = GaussTri(order);
      R{1} = GaussTet(order);
    end
    function R = isFeasible(points)
      tol = 1e-12;
      R = (all(points>-tol, 2) & 1-sum(points,2)>-tol);
    end
    function R = getCenterLoc()
      R = [1 1 1]/4;
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