classdef Algorithm < SOFE
  properties
    pde
  end
  methods % constructor
    function obj = Algorithm(pde)
      obj.pde = pde;
    end
  end
end
