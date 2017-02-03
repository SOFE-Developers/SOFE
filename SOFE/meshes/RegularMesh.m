classdef RegularMesh < Mesh
  methods % constructor
    function obj = RegularMesh(NVec, diam, isTri)
      fprintf('Generate mesh ... ');
      dim = numel(NVec);
      GRID = cell(1,dim);
      for d = 1:dim
        GRID{d} = linspace(diam(d,1),diam(d,2),NVec(d)+1);
      end
      [nodes, elems] = Mesh.getTensorProductMesh(GRID, isTri);
      obj = obj@Mesh(nodes, elems);
      fprintf('DONE\n');
    end
  end
end