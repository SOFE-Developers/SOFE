classdef IterativeSolver < Solver
  properties
    type
    precon
    tol = 1e-10;
    maxit = 1000;
  end
  methods % constructor
    function obj = IterativeSolver(pde, type, precon)
      obj = obj@Solver(pde);
      obj.type = type;
      obj.precon = precon;
    end
  end
  methods % solve
    function R = solve(obj, A, b, varargin) % [x0]
      switch obj.precon
        case 'diag'
          M1 = spdiags(diag(A), 0, size(A,1), size(A,1));
          M2 = speye(size(M1)); %#ok<NASGU>
        case 'ilu'
          [M1,M2] = ilu(A); %#ok<ASGLU>
        case 'ichol'
          M1 = ichol(A,struct('michol','on'));
          M2 = M1'; %#ok<NASGU>
        case 'none'
          M1 = speye(size(A,1)); M2 = M1; %#ok<NASGU>
        otherwise
          warning('Unknown preconditioner');
      end
      try
        eval(['R = ' obj.type '(A, b, obj.tol, obj.maxit, M1, M2, varargin{:});']);
      catch
        warning('Solving failed!');
        R = zeros(size(b));
      end
    end
    function R = solve2(obj, b, freeI, freeJ)
      afun = @(x)obj.pde.applySystem(x, freeI, freeJ);
      R = bicgstab(afun, b, obj.tol, 74);
      R = R(freeI);
    end
  end
end
