classdef VisualizerMesh1D < VisualizerMesh
  methods % constructor
    function obj = VisualizerMesh1D(mesh)
      obj = obj@VisualizerMesh(mesh);
    end
  end
  methods % display
    function show(obj)
      switch size(obj.mesh.nodes, 2)
        case 1
          plot(obj.mesh.nodes, zeros(obj.mesh.getNumber(0),1), '*');
        case 2
          plot(obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), '*');
        case 3
          plot3(obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), obj.mesh.nodes(:,3), '*');
      end
    end
    function showEntity(obj, dim)
      center = obj.mesh.getCenter(dim);
      nE = obj.mesh.getNumber(dim);
      switch dim
        case 1
          color = [0 1 0];
        case 0
          color = [1 0 1];
      end
      text(center, 0*center, num2str((1:nE)'),'Color',color,'FontSize', 18);
    end
  end
end