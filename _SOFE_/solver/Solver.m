classdef Solver < SOFEClass
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
    function solve(obj)
    end
  end
end
