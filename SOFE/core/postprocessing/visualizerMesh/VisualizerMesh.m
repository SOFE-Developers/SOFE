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
      R = [];
      dim = mesh.dimW;
      eval(['R = VisualizerMesh' num2str(dim) 'D(mesh);']);
    end
  end
end