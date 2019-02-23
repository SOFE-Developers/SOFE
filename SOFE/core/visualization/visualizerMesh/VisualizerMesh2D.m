classdef VisualizerMesh2D < VisualizerMesh
  methods % constructor
    function obj = VisualizerMesh2D(mesh)
      obj = obj@VisualizerMesh(mesh);
    end
  end
  methods % display
    function show(obj)
      elem = obj.mesh.topology.getEntity('0');
      if obj.mesh.topology.isSimplex, I = [1 2 3]; else, I = [1 2 4 3]; end
      if size(obj.mesh.nodes, 2) == 2
        h = trisurf(elem(:,I), obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), zeros(size(obj.mesh.nodes,1),1));
        view(0,90);
      else
        h = trimesh(elem(:,I), obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), obj.mesh.nodes(:,3));
      end
      set(h,'facecolor',[0.5 0.8 0.5],'edgecolor','k','LineWidth',1.5);
      axis equal
    end
    function showEntity(obj, dim)
      center = obj.mesh.getCenter(dim);
      nE = obj.mesh.topology.getNumber(dim);
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
end