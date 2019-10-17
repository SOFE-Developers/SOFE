classdef LsqMinSol < Solver
  methods % constructor
    function obj = LsqMinSol()
      obj = obj@Solver();
    end
  end
  methods % solve
    function R = solve(obj, A, b, I, J, shift)
      assert(obj.pde.createSys, 'System must be created');
      b = b - A*shift;
      b = b(I);
      A = A(I,J);
      R = zeros(size(shift));
      R(J) = lsqminnorm(A, b);
      R = shift + R;
    end
  end
end
