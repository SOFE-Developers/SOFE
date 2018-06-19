classdef Algorithm < SOFE
  properties
    pde
    solver
  end
  methods % constructor
    function obj = Algorithm(pde, varargin) % [solver]
      obj.pde = pde;
      if ~isempty(varargin)
        obj.setSolver(varargin{1});
      end
    end
    function setSolver(obj, solver)
      obj.solver = solver;
      obj.solver.setPDE(obj.pde);
    end
  end
end
