classdef IterativeSolver < Solver
  properties
    type
    precon
    tol = 1e-10;
    maxit = 1000;
  end
  methods % constructor
    function obj = IterativeSolver(pde, type, varargin) % [precon]
      obj = obj@Solver(pde);
      obj.type = type;
      if ~isempty(varargin)
        obj.precon = varargin{1};
      else
        obj.precon = 'none';
      end
    end
  end
  methods % solve
    function R = compute(obj)
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      [freeI, freeJ] = obj.pde.getFreeDoFs();
      fprintf('%d DoFs\n', sum(freeJ));
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);    
      t = tic; obj.output('Begin solve ...', 1);
      M1 = []; M2 = [];
      if obj.pde.createSys 
        b = obj.pde.loadVec - obj.pde.stiffMat*obj.pde.getShift();
        A = obj.pde.stiffMat(freeI, freeJ);
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
        b = obj.pde.loadVec - obj.pde.applySystem(obj.pde.getShift()); 
        A = @(x)obj.pde.applySystem(x, freeI, freeJ);
      end
      b = b(freeI); %#ok<*NASGU>
      R = obj.pde.getShift();
      R(freeJ) = eval([obj.type '(A, b, obj.tol, obj.maxit, M1, M2);']);
      obj.solution = R;
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
