classdef FixPoint < Algorithm
  properties
    nIt
    history
  end
  methods % constructor
    function obj = FixPoint(pde, nIt, varargin) % [solver]
      obj = obj@Algorithm(pde, varargin{:});
      obj.nIt = nIt;
      obj.history = cell(nIt,1);
    end
  end
  methods % integrate
    function compute(obj)
      for k = 1:obj.nIt
        obj.pde.assemble();
        [freeI, freeJ] = obj.pde.getFreeDoFs();
        shift = obj.pde.getShift();
        obj.history{k} = obj.solver.solve(obj.pde.stiffMat, obj.pde.loadVec, freeI, freeJ, shift);
        diff = norm(obj.history{k} - obj.pde.state);
        obj.pde.setState(0.0, obj.history{k});
        fprintf('step: %d / %d, diff: %d\n', k, obj.nIt, diff);
      end
    end
  end
end