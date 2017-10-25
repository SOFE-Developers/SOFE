classdef VisualizerMeshInt < VisualizerMesh
  methods % constructor
    function obj = VisualizerMeshInt(mesh)
      obj = obj@VisualizerMesh(mesh);
    end
  end
  methods % display
    function show(obj)
      top = obj.mesh.topology;
      switch size(top.nodes, 2)
        case 1
          plot(top.nodes, zeros(top.getNumber(0),1), '*');
        case 2
          plot(top.nodes(:,1), top.nodes(:,2), '*');
        case 3
          plot3(top.nodes(:,1), top.nodes(:,2), top.nodes(:,3), '*');
      end
    end
    function showEntity(obj, dim)
      top = obj.mesh.topology;
      center = top.getCenter(dim);
      nE = top.getNumber(dim);
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