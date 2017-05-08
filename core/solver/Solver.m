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
    function reset(obj) %#ok<MANU>
    end
    function solve(obj, A, b) %#ok<INUSD>
    end
  end
end
