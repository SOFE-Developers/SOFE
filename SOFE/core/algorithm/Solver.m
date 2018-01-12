classdef Solver < Algorithm
  properties
    solution
  end
  methods % constructor
    function obj = Solver(pde)
      obj = obj@Algorithm(pde);
    end
  end
  methods
    function R = getSolution(obj, idx)
      R = obj.solution(obj.pde.J(idx,1):obj.pde.J(idx,2));
    end
    function setSolution(obj, sol, idx)
      if size(obj.solution,1)<obj.pde.J(idx,2)
        obj.solution = [obj.solution;zeros(diff(obj.pde.J(idx,:))+1,1)];
      end
      obj.solution(obj.pde.J(idx,1):obj.pde.J(idx,2)) = sol;
    end
  end
end
