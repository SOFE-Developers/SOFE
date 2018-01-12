classdef FixPointSolver < Solver
  properties
    nIt
    linSolver
  end
  methods % constructor
    function obj = FixPointSolver(p, nIt, varargin) % [initCond]
      obj = obj@Solver(p);
      obj.nIt = nIt;
      obj.linSolver = DirectSolver(p);
      if ~isempty(varargin)
        obj.pde.setState(varargin{1});
      end
    end
  end
  methods % integrate
    function compute(obj)
      obj.solution = zeros(obj.pde.nDoF,1);
      N = obj.pde.fesTrial{1}.getNDoF();
      for k = 1:obj.nIt
        obj.linSolver.compute();
        obj.pde.setState(0.0, obj.linSolver.solution);
        diff = norm(obj.solution(1:N) - obj.linSolver.solution(1:N));
        fprintf('step: %d / %d, diff: %d\n', k, obj.nIt, diff);
        obj.solution = obj.linSolver.solution;
      end
    end
  end
end