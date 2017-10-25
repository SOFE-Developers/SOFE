classdef VisualizerMeshTri < VisualizerMesh
  methods % constructor
    function obj = VisualizerMeshTri(mesh)
      obj = obj@VisualizerMesh(mesh);
    end
  end
  methods % display
    function show(obj)
      top = obj.mesh.topology;
      elem = top.getEntity('0');
      if size(top.nodes, 2) == 2
        h = trisurf(elem, top.nodes(:,1), top.nodes(:,2), zeros(size(top.nodes,1),1));
        view(0,90);
      else
        h = trimesh(elem, top.nodes(:,1), top.nodes(:,2), top.nodes(:,3));
      end
      set(h,'facecolor',[0.5 0.8 0.5], 'edgecolor', 'k');
    end
    function showEntity(obj, dim)
      top = obj.mesh.topology;
      center = obj.mesh.getCenter(dim);
      nE = top.getNumber(dim);
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