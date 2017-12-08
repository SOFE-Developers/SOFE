classdef DirectSolver < Solver
  properties
    isLU = false;
    matDecomp
  end
  methods % constructor
    function obj = DirectSolver(pde, varargin) % [isLU]
      obj = obj@Solver(pde);
      if nargin > 1, obj.isLU = varargin{1}; end
    end
    function reset(obj)
      obj.matDecomp = [];
    end
  end
  methods % solve
    function R = solve(obj, A, b, varargin)
      if obj.isLU
        if isempty(obj.matDecomp)
          [obj.matDecomp.L,obj.matDecomp.U,obj.matDecomp.P,obj.matDecomp.Q,obj.matDecomp.R] = lu(A);
        end
        R = obj.matDecomp.Q*(obj.matDecomp.U\(obj.matDecomp.L\(obj.matDecomp.P*(obj.matDecomp.R\b))));
      else
        R = A\b;
      end
    end
  end
end
