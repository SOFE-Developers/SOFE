classdef Solver < SOFE
  properties
    pde
  end
  methods % constructor
    function obj = Solver(pde)
      obj.pde = pde;
    end
  end
  methods
    function reset(obj)
    end
    function solve(obj, A, b)
    end
  end
end
