classdef IterativeSolver2 < Solver
  properties
    type
    precon
    tol = 1e-10;
    maxit = 1000;
  end
  methods % constructor
    function obj = IterativeSolver2(pde, type, precon)
      obj = obj@Solver(pde);
      obj.type = type;
      obj.precon = precon;
    end
  end
  methods % solve
    function R = compute(obj)
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      fprintf('%d DoFs\n', sum(obj.pde.fDoFsTrial));
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);    
      t = tic; obj.output('Begin solve ...', 1);
      M1 = []; M2 = [];
      if obj.pde.createSys 
        b = obj.pde.loadVec - obj.pde.stiffMat*obj.pde.shift;
        A = obj.pde.stiffMat(obj.pde.fDoFsTest, obj.pde.fDoFsTrial);
        switch obj.precon
          case 'diag'
            M1 = spdiags(diag(A), 0, size(A,1), size(A,1));
            M2 = speye(size(M1));
          case 'ilu'
            [M1,M2] = ilu(A); %#ok<ASGLU>
          case 'ichol'
            M1 = ichol(A,struct('michol','on'));
            M2 = M1';
          case 'none'
          otherwise
            warning('Unknown preconditioner');
        end
      else
        b = obj.pde.loadVec - obj.pde.applySystem(obj.pde.shift); 
        A = @(x)obj.pde.applySystem(x, 1);
      end
      b = b(obj.pde.fDoFsTest); %#ok<*NASGU>
      R = obj.pde.shift;
      R(obj.pde.fDoFsTrial) = eval([obj.type '(A, b, obj.tol, obj.maxit, M1, M2);']);
      obj.pde.solution = R;
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
