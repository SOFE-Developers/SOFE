classdef MeshTopologyInt < MeshTopology
  methods % constructor
    function obj = MeshTopologyInt(elem)
      obj = obj@MeshTopology(1);
      obj.update(elem);
    end
    function update(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{1,1} = (1:max(obj.connectivity{2,1}(:)))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
    end
  end
  methods % refinement
    function P = uniformRefine(obj)
      el = obj.getEntity(1);
      nE = obj.getNumber(1); nN = obj.getNumber(0);     
      P = [speye(nN); sparse(repmat((1:nE)',1,2), el, 0.5)];
      newIndices = nN + (1:nE);
      el = [el newIndices'];
      el = [el(:,[1 3]); el(:,[3 2])];
      obj.update(el);
    end
  end
  methods(Static = true)
    function R = getQuadRule(order)
      R{2} = GaussPoint();
      R{1} = GaussInt(order);
    end
    function R = isFeasible(points, tol)
      R = (points>-tol & points<1+tol);
    end
    function R = getCenterLoc()
      R = 1/2;
    end
  end
end