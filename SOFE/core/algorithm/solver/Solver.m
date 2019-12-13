classdef Solver < SOFE
  properties
    pde
  end
  methods % constructor
    function obj = Solver()
      obj = obj@SOFE();
    end
    function setPDE(obj, pde)
      obj.pde = pde;  
    end
  end
end
