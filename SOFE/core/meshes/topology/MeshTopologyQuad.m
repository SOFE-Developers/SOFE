classdef MeshTopologyQuad < MeshTopology
  methods % constructor
    function obj = MeshTopologyQuad(elem, dimP)
      obj = obj@MeshTopology(dimP);
      obj.update(elem);
    end
    function update(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{2,1} = [elem(:,[1,2]); elem(:,[3,4]); elem(:,[1,3]); elem(:,[2,4])];
      [obj.connectivity{2,1}, ~, e2F] = unique(sort(obj.connectivity{2,1},2),'rows'); 
      obj.connectivity{3,2} = reshape(e2F, size(elem,1), []);
      %
      obj.connectivity{1,1} = (1:max(obj.connectivity{2,1}(:)))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
    end
  end
  methods % connectivity information   
    function R = getOrientation(obj, varargin)
      e = obj.getEntity('0');
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
    function P = uniformRefine(obj)
      faces = obj.getEntity(1);
      el = obj.getEntity(2);
      nE = obj.getNumber(2); nF = obj.getNumber(1); nN = obj.getNumber(0);
      P = [eye(nN); sparse(repmat((1:nF)',1,2), faces, 0.5); ...
                    sparse(repmat((1:nE)',1,4), el, 0.25)];
      newIndicesF = nN + (1:nF);
      newIndicesE = nN + nF + (1:nE)';
      el = [el newIndicesF(obj.connectivity{3,2}) newIndicesE];
      el = [el(:,[1 5 7 9]);el(:,[5 2 9 8]);el(:,[7 9 3 6]);el(:,[9 8 6 4])];
      obj.update(el);
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