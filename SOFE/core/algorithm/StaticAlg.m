classdef StaticAlg < Algorithm
  properties
    solution
  end
  methods % constructor
    function obj = StaticAlg(pde, varargin) % [solver]
      obj = obj@Algorithm(pde, varargin{:});
    end
  end
  methods
    function R = compute(obj)
      assert(~isempty(obj.solver), 'Solver not set!');
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      [freeI, freeJ] = obj.pde.getFreeDoFs();
      shift = obj.pde.getShift();
      fprintf('%d DoFs\n', size(obj.pde.loadVec,1));
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);
      %
      t = tic; obj.output('Begin solve ...', 1);
      obj.solution = obj.solver.solve(obj.pde.stiffMat, obj.pde.loadVec, freeI, freeJ, shift);
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
