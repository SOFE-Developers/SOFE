classdef MeshTopologyInt < MeshTopology
  methods % constructor
    function obj = MeshTopologyInt(nodes, elem, dimP)
      obj = obj@MeshTopology(nodes, dimP);
      obj.updateConnectivity(elem);
    end
    function updateConnectivity(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
    end
  end
  methods % refinement
    function P = uniformRefine(obj)
      el = obj.getEntity(1);
      nE = obj.getNumber(1); nN = obj.getNumber(0);     
      P = [eye(nN); sparse(repmat((1:nE)',1,2), el, 0.5)];
      obj.nodes = P*obj.nodes;
      newIndices = nN + (1:nE);
      el = [el newIndices'];
      el = [el(:,[1 3]); el(:,[3 2])];
      obj.updateConnectivity(el);
    end
  end
  methods(Static = true)
    function R = getQuadRule(order)
      R{2} = GaussPoint();
      R{1} = GaussInt(order);
    end
    function R = isFeasible(points)
      tol = 1e-12;
      R = (points>-tol & points<1+tol);
    end
    function R = getCenterLoc()
      R = 1/2;
    end
  end
end