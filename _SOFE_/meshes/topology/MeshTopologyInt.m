classdef MeshTopologyInt < MeshTopology
  methods % constructor
    function obj = MeshTopologyInt(nodes, elem, dimP)
      obj = obj@MeshTopology(nodes, elem, dimP);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
    end
  end
  methods % display
    function show(obj)
      switch size(obj.nodes, 2)
        case 1
          plot(obj.nodes, zeros(obj.getNumber(0),1), '.');
        case 2
          plot(obj.nodes(:,1), obj.nodes(:,2), '.');
        case 3
          plot3(obj.nodes(:,1), obj.nodes(:,2), obj.nodes(:,3), '.');
      end
    end
    function showNodeVector(obj, U)
      plot(obj.nodes, U(1:obj.getNumber(0)));
    end
  end
  methods % refinement
    function uniformRefine(obj)
      % data
      el = obj.getEntity(1);
      nE = obj.getNumber(1); nN = obj.getNumber(0);
      % node coords
      obj.nodes = [obj.nodes; obj.getCenter(0)];
      % node indices
      newIndices = nN + (1:nE);
      el = [el newIndices'];
      % elem
      obj.connectivity{2,1} = [el(:,[1 3]); el(:,[3 2])];
      % update
      obj.updateConnectivity();
      %
      obj.notifyObservers();
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
  end
end