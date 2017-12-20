classdef DirectSolver2 < Solver
  properties
    isCachingLU = false;
    matDecomp
  end
  methods % constructor
    function obj = DirectSolver2(pde, varargin) % [isLU]
      obj = obj@Solver(pde);
      if nargin > 1, obj.isCachingLU = varargin{1}; end
    end
    function notify(obj)
      obj.matDecomp = [];
    end
  end
  methods % solve
    function R = compute(obj)
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      fprintf('%d DoFs\n', sum(obj.pde.fDoFsTrial));
      %
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);    
      t = tic; obj.output('Begin solve ...', 1);
      if ~obj.pde.createSys 
        error('System must be created');
      end
      b = obj.pde.loadVec - obj.pde.stiffMat*obj.pde.shift;
      b = b(obj.pde.fDoFsTest);
      A = obj.pde.stiffMat(obj.pde.fDoFsTest, obj.pde.fDoFsTrial);
      R = obj.pde.shift;
      if obj.isCachingLU
        if isempty(obj.matDecomp)
          [obj.matDecomp.L,obj.matDecomp.U,obj.matDecomp.P,obj.matDecomp.Q,obj.matDecomp.R] = lu(A);
        end
        R(obj.pde.fDoFsTrial) = obj.matDecomp.Q*(obj.matDecomp.U\(obj.matDecomp.L\(obj.matDecomp.P*(obj.matDecomp.R\b))));
      else
        R(obj.pde.fDoFsTrial) = A\b;
      end
      obj.pde.solution = R;
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
