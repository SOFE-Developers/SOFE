classdef Solver < SOFE
  properties
    pde
  end
  methods % constructor
    function obj = Solver(pde)
      obj.pde = pde;
    end
  end
end
