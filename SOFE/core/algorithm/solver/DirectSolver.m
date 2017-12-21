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
      fprintf('%d DoFs\n', sum(obj.pde.fTrial));
      %
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);    
      t = tic; obj.output('Begin solve ...', 1);
      if ~obj.pde.createSys 
        error('System must be created');
      end
      b = obj.pde.loadVec - obj.pde.stiffMat*obj.pde.shift;
      b = b(obj.pde.fTest);
      A = obj.pde.stiffMat(obj.pde.fTest, obj.pde.fTrial);
      R = obj.pde.shift;
      if obj.isCaching
        if isempty(obj.cache)
          [obj.cache.L,obj.cache.U,obj.cache.P,obj.cache.Q,obj.cache.R] = lu(A);
        end
        R(obj.pde.fTrial) = obj.cache.Q*(obj.cache.U\(obj.cache.L\(obj.cache.P*(obj.cache.R\b))));
      else
        R(obj.pde.fTrial) = A\b;
      end
      obj.pde.solution = R;
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
