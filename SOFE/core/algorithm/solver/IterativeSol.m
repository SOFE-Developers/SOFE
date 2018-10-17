classdef IterativeSol < Solver
  properties
    type
    precon
    tol = 1e-10;
    maxit = 1000;
  end
  methods % constructor
    function obj = IterativeSol(type, precon)
      obj = obj@Solver();
      obj.type = type;
      obj.precon = precon;
    end
  end
  methods % solve
    function R = solve(obj, A, b, I, J, shift)
      M1 = []; M2 = [];
      b = b - A*shift;
      A = A(I, J);
      switch obj.precon
        case 'diag'
          M1 = spdiags(diag(A), 0, size(A,1), size(A,1));
          M2 = speye(size(M1));
        case 'ilu'
          [M1,M2] = ilu(A); %#ok<ASGLU>
        case 'ichol'
          M1 = ichol(A,struct('michol','on'));
          M2 = M1';
        case {[], 'none', 'no'}
        otherwise
          warning('Unknown preconditioner');
      end
      b = b(I); %#ok<*NASGU>
      R = zeros(size(shift));
      switch obj.type
        case 'gmres'
          nR = 40;
          R(J) = eval([obj.type '(A, b,' num2str(nR) ', obj.tol, obj.maxit, M1, M2);']);
        otherwise
          R(J) = eval([obj.type '(A, b, obj.tol, obj.maxit, M1, M2);']);
      end
      R = shift + R;
    end
  end
end
