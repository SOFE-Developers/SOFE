classdef VisualizerMesh < SOFE
  properties
    mesh
  end
  methods % constructor & more
    function obj = VisualizerMesh(mesh)
      obj.mesh = mesh;
    end
  end
  methods(Static = true)
    function R = create(mesh)
      name = class(mesh.topology);
      name = name(13:end);
      eval(['R = VisualizerMesh' name '(mesh);']);
    end
  end
end