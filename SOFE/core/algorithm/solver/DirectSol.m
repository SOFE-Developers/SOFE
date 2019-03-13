classdef DirectSol < Solver
  properties
    isLU
    cache
  end
  methods % constructor
    function obj = DirectSol(varargin) % [isLU]
      obj = obj@Solver();
      if ~isempty(varargin), obj.isLU = varargin{1}; else, obj.isLU = 0; end
    end
  end
  methods % solve
    function R = solve(obj, A, b, I, J, shift)
      assert(obj.pde.createSys, 'System must be created');
      b = b - A*shift;
      b = b(I);
      A = A(I,J);
      R = zeros(size(shift));
      if obj.isLU
        if isempty(obj.cache)
          [obj.cache.L,obj.cache.U,obj.cache.P,obj.cache.Q,obj.cache.R] = lu(A);
        end
        R(J) = obj.cache.Q*(obj.cache.U\(obj.cache.L\(obj.cache.P*(obj.cache.R\b))));
      else
        R(J) = A\b;
      end
      R = shift + R;
    end
  end
end
