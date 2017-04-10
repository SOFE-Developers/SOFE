classdef MeshTopologyTri < MeshTopology
  methods % constructor
    function obj = MeshTopologyTri(nodes, elem, dimP)
      elem = MeshTopologyTri.renumber(nodes, elem);
      obj = obj@MeshTopology(nodes, elem, dimP);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectivity{2,1} = [obj.connectivity{3,1}(:,[1,2]); ...
                               obj.connectivity{3,1}(:,[2,3]); ...
                               obj.connectivity{3,1}(:,[1,3])];
      [obj.connectivity{2,1}, dmy, e2F] = unique(sort(obj.connectivity{2,1},2),'rows');      
      obj.connectivity{3,2} = reshape(e2F, [], 3);
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
      R(e(:,2)>e(:,3),2) = -1;
      R(e(:,1)>e(:,3),3) = -1;
      if nargin > 3
        R = R(varargin{3},:);
      end
    end
    function R = getNormalOrientation(obj, varargin)
      R = obj.getOrientation();
      R(:,3) = -R(:,3);
    end
    function R = upliftPoints(obj, points, fLoc, orient)
      zz = zeros(size(points));
      if orient<0
        points = 1-points;
      end
      switch fLoc
        case 1
          R = [points, zz];
        case 2
          R = [1-points, points];
        case 3
          R = [zz, points];
      end
    end
  end
  methods % mesh information
    function R = getMeasure(obj, dim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      ee = obj.getEntity(dim); ee = ee(I,:);
      nodes = obj.getEntity(0);
      switch dim
        case 2 % element
          v1 = nodes(ee(:,2),:) - nodes(ee(:,1),:);
          v2 = nodes(ee(:,3),:) - nodes(ee(:,1),:);
          R = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))/2;
        case 1 % face
          v = nodes(ee(:,2),:) - nodes(ee(:,1),:);
          R = sum(v.^2,2).^0.5;
       end
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(2);
      nE = obj.getNumber(2); nF = obj.getNumber(1); nN = obj.getNumber(0);
      % node coords
      obj.nodes = [obj.nodes; obj.getCenter(1)];
      % node indices
      newIndices = nN + (1:nF);
      el = [el newIndices(obj.connectivity{3,2})];
      % elem
      obj.connectivity{3,1} = [el(:,[1 4 6]);el(:,[4 2 5]);el(:,[6 5 3]);el(:,[5 6 4])];
      % update
      obj.updateConnectivity();
      %
      obj.notifyObservers();
    end
  end
  methods % display
    function show(obj)
      elem = obj.getEntity(obj.dimP);
      nodes = obj.getEntity(0);
      if size(nodes, 2) == 2
        trisurf(elem, nodes(:,1), nodes(:,2), zeros(size(nodes,1),1));
        view(0,90);
      else
        h = trimesh(elem, nodes(:,1), nodes(:,2), nodes(:,3));
        set(h, 'EdgeColor','black')
      end
    end
    function showEntity(obj, dim) % for debug
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
    function R = renumber(nodes, elem)
      if isempty(elem), R = []; return; end
      v1 = nodes(elem(:,2),:) - nodes(elem(:,1),:);
      v2 = nodes(elem(:,3),:) - nodes(elem(:,1),:);
      I = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))/2;
      I = I<0;
      if any(I)
        fprintf('Elements renumbered!\n');
        elem(I,:) = elem(I, [2 1 3]);
      end
      R = elem;
    end
    function R = isFeasible(points)
      tol = 1e-12;
      R = (all(points>-tol, 2) & 1-sum(points,2)>-tol);
    end
    function R = getQuadRule(quadOrder)
      R{3} = GaussPoint();
      R{2} = GaussInt(quadOrder);
      R{1} = GaussTri(quadOrder);
    end
    function R = getBarycenterRef()
      R = [1 1]/3;
    end
    function R = getFacePointList(pointsFace)
      nO = 2;
      R = cell(nO,1);
      dualPoints = 1 - pointsFace;
      zeroPoints = zeros(size(pointsFace));
      % face 1
      R{1}(1,:,:) = [pointsFace zeroPoints];
      R{2}(1,:,:) = [dualPoints zeroPoints];
      % face 2
      R{1}(2,:,:) = [dualPoints pointsFace];
      R{2}(2,:,:) = [pointsFace dualPoints];
      % face 3
      R{1}(3,:,:) = [zeroPoints pointsFace];
      R{2}(3,:,:) = [zeroPoints dualPoints];
    end
    function R = getEquiPoints(order)
      ls = linspace(0,1, order+1);
      [X, Y] = meshgrid(ls, ls);
      R = [X(:) Y(:)];
      R = R(sum(R,2)<=1,:);
    end
    function R = getRefPointsChild(P, childNr)
      switch childNr
        case 0
          R = P;
        case 1
          R = 0.5*P;
        case 2
          R = 0.5*P; R(:,1) = R(:,1) + 0.5;
        case 3
          R = 0.5*P; R(:,2) = R(:,2) + 0.5;
        case 4
          R = -0.5*P; R = R + 0.5;
      end
    end
  end
end