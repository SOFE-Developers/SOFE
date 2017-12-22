classdef DirectSolver < Solver
  properties
    isCaching = false;
    cache
  end
  methods % constructor
    function obj = DirectSolver(pde, varargin) % [isLU]
      obj = obj@Solver(pde);
      if nargin > 1, obj.isCaching = varargin{1}; end
    end
    function notify(obj)
      obj.cache = [];
    end
  end
  methods % solve
    function R = compute(obj)
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      [freeI, freeJ] = obj.pde.getFreeDoFs();
      fprintf('%d DoFs\n', sum(freeJ));
      %
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);    
      t = tic; obj.output('Begin solve ...', 1);
      if ~obj.pde.createSys 
        error('System must be created');
      end
      b = obj.pde.loadVec - obj.pde.stiffMat*obj.pde.getShift();
      b = b(freeI);
      A = obj.pde.stiffMat(freeI, freeJ);
      R = obj.pde.getShift();
      if obj.isCaching
        if isempty(obj.cache)
          [obj.cache.L,obj.cache.U,obj.cache.P,obj.cache.Q,obj.cache.R] = lu(A);
        end
        R(freeJ) = obj.cache.Q*(obj.cache.U\(obj.cache.L\(obj.cache.P*(obj.cache.R\b))));
      else
        R(freeJ) = A\b;
      end
      obj.solution = R;
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
