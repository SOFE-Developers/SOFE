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
      [obj.connectivity{2,1}, ~, e2F] = unique(sort(obj.connectivity{2,1},2),'rows'); 
      obj.connectivity{3,2} = reshape(e2F, [], 4);
      %
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
    end
  end
  methods % connectivity information   
    function R = getOrientation(obj, varargin)
      e = obj.getEntity(obj.dimP);
      R = ones(size(e));
      R(e(:,1)>e(:,2),1) = -1;
      R(e(:,3)>e(:,4),2) = -1;
      R(e(:,1)>e(:,3),3) = -1;
      R(e(:,2)>e(:,4),4) = -1;
      if nargin > 3
        R = R(varargin{3},:);
      end
    end
    function R = getNormalOrientation(obj, varargin)
      R = obj.getOrientation();
      R(:,[2 3]) = -R(:,[2 3]);
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(2);
      nE = obj.getNumber(2); nF = obj.getNumber(1); nN = obj.getNumber(0);
      obj.nodes = [obj.nodes; obj.getCenter(1); obj.getCenter(2)];
      newIndicesF = nN + (1:nF);
      newIndicesE = nN + nF + (1:nE)';
      el = [el newIndicesF(obj.connectivity{3,2}) newIndicesE];
      obj.connectivity{3,1} = [el(:,[1 5 7 9]);el(:,[5 2 9 8]);el(:,[7 9 3 6]);el(:,[9 8 6 4])];
      obj.updateConnectivity();
      obj.notifyObservers();
    end
  end
  methods % display
    function show(obj)
      elem = obj.getEntity(obj.dimP);
      if size(obj.nodes, 2) == 2
        trisurf(elem(:,[1 2 4 3]), obj.nodes(:,1), obj.nodes(:,2), zeros(size(obj.nodes,1),1));
        view(0,90);
      else
        h = trimesh(elem(:,[1 2 4 3]), obj.nodes(:,1), obj.nodes(:,2), obj.nodes(:,3));
        set(h, 'EdgeColor','black')
      end
    end
    function showEntity(obj, dim)
      center = obj.getCenter(dim);
      nE = obj.getNumber(dim);
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
  end
  methods(Static = true)
    function R = getQuadRule(order)
      R{3} = GaussPoint();
      R{2} = GaussInt(order);
      R{1} = GaussQuad(order);
    end
    function R = isFeasible(points)
      tol = 1e-12;
      R = (points(:,1)>-tol & points(:,1)<1+tol & points(:,2)>-tol & points(:,2)<1+tol);
    end
    function R = getCenterLoc()
      R = [1 1]/2;
    end
    function R = upliftPoints(points, fLoc, orient)
      zz = zeros(size(points)); oo = ones(size(points));
      if orient<0
        points = 1-points;
      end
      switch fLoc
        case 1
          R = [points, zz];
        case 2
          R = [points, oo];
        case 3
          R = [zz, points];
        case 4
          R = [oo, points];
      end
    end
  end
end