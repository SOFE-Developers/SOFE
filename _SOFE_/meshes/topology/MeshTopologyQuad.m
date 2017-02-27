classdef MeshTopologyQuad < MeshTopology
  methods % constructor
    function obj = MeshTopologyQuad(nodes, elems, dimP)
      obj = obj@MeshTopology(nodes, elems, dimP);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectivity{2,1} = [obj.connectivity{3,1}(:,[1,2]); ...
                               obj.connectivity{3,1}(:,[3,4]); ...
                               obj.connectivity{3,1}(:,[1,3]); ...
                               obj.connectivity{3,1}(:,[2,4])];
      [obj.connectivity{2,1}, dmy, e2F] = unique(sort(obj.connectivity{2,1},2),'rows'); 
      obj.connectivity{3,2} = reshape(e2F, [], 4);
      %
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
    end
  end
  methods % connectivity information   
    function R = getElem2Face(obj)
      R = obj.connectivity{3,2}; 
    end
    function R = getFace2Elem(obj)
      orient = obj.getOrientation();
      nE = obj.getNumber(2);
      nF = obj.getNumber(1);
      col = repmat([1 2 2 1], nE, 1);
      col(orient<0) = -col(orient<0)+3;
      R = full(sparse(obj.getElem2Face(), col, repmat((1:nE)',1,4),nF,2));
    end
    function R = getFaceType(obj)
      orient = obj.getOrientation();
      nE = obj.getNumber(2);
      nF = obj.getNumber(1);
      col = repmat([1 2 2 1], nE, 1);
      col(orient<0) = -col(orient<0)+3;
      R = full(sparse(obj.getElem2Face(), col, ones(nE,1)*[1 2 3 4],nF,2));
    end
    function R = getOrientation(obj, varargin)
      e = obj.getEntity(2);
      R = ones(size(e));
      R(e(:,1)>e(:,2),1) = -1;
      R(e(:,3)>e(:,4),2) = -1;
      R(e(:,1)>e(:,3),3) = -1;
      R(e(:,2)>e(:,4),4) = -1;
    end
    function R = getNormalOrientation(obj, varargin)
      R = obj.getOrientation();
      R(:,[2 3]) = -R(:,[2 3]);
    end
  end
  methods % mesh information
    function R = getQuadRule(obj, order)
      R{3} = GaussPoint();
      R{2} = GaussInt(order);
      R{1} = GaussQuad(order);
    end
    function R = getBarycenterRef(obj)
      R = [1 1]/2;
    end
    function R = getMeasure(obj, dim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      ee = obj.getEntity(dim); ee = ee(I,:);
      switch dim
        case 2 % element
          v1 = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          v2 = obj.nodes(ee(:,3),:) - obj.nodes(ee(:,1),:);
          R = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1));
          if all(all(obj.nodes(ee(:,1),:) + v1+v2 - obj.nodes(ee(:,4),:)))
            warning('! Area only valid for parallelograms !');
          end
        case 1 % face
          v = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          R = sum(v.^2,2).^0.5;
       end
    end
    function R = isFeasible(obj, points)
      tol = 1e-12;
      R = (points(:,1)>-tol & points(:,1)<1+tol & points(:,2)>-tol & points(:,2)<1+tol);
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(2);
      nE = obj.getNumber(2); nF = obj.getNumber(1); nN = obj.getNumber(0);
      % node coords
      obj.nodes = [obj.nodes; obj.getCenter(1); obj.getCenter(2)];
      % node indices
      newIndicesF = nN + (1:nF);
      newIndicesE = nN + nF + (1:nE)';
      el = [el newIndicesF(obj.connectivity{3,2}) newIndicesE];
      % elem
      obj.connectivity{3,1} = [el(:,[1 5 7 9]);el(:,[5 2 9 8]);el(:,[7 9 3 6]);el(:,[9 8 6 4])];
      % update
      obj.updateConnectivity();
      %
      obj.notifyObservers();
    end
  end
  methods % display
    function show(obj)
      elem = obj.getEntity(2);
      if size(obj.nodes, 2) == 2
        trisurf(elem(:,[1 2 4 3]), obj.nodes(:,1), obj.nodes(:,2), zeros(size(obj.nodes,1),1));
        view(0,90);
      else
        h = trimesh(elem(:,[1 2 4 3]), obj.nodes(:,1), obj.nodes(:,2), obj.nodes(:,3));
        set(h, 'EdgeColor','black')
      end
    end
    function showEntity(obj, dim) % for debug
      center = obj.getCenter(dim);
      if dim == 0
        nE = size(obj.nodes, 1);
      else
        nE = obj.getNumber(dim);
      end
      switch dim
        case 2
          color = [1 0 1];
        case 1
          color = [1 1 0];
        case 0
          color = [0 0.5 1];
      end
      text(center(:,1), center(:,2), num2str((1:nE)'),'Color',color,'FontSize', 18);
    end
    function showNodeVector(obj, U)
      elem = obj.getEntity(2);
      h = trimesh(elem(:,[1 2 4 3]), obj.nodes(:,1), obj.nodes(:,2), U(1:obj.getNumber(0)));
      set(h,'facecolor','interp','edgecolor','k');
    end
  end
end