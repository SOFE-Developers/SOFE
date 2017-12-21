classdef Solver < Algorithm
  properties
    solution
  end
  methods % constructor
    function obj = Solver(pde)
      obj = obj@Algorithm(pde);
    end
  end
end
