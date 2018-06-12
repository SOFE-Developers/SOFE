classdef MultiStepMethod < TimeStep
  methods % constructor
    function obj = MultiStepMethod(data, M0, pde, solver)
      obj = obj@TimeStep(data, M0, pde, solver);
      assert(obj.pde.nArgIn.coeff<2, 'Data (still) must be independent of time');
      obj.nS = max(numel(obj.data.alpha), numel(obj.data.beta))-1;
      obj.nK = 1;
    end
  end
  methods % computation
    function R = compute(obj, I, u0)
      obj.pde.setState(I(2), u0);
      obj.pde.assemble();
      lhs = obj.data.alpha(1)*obj.M0.stiffMat + (diff(I)*obj.data.beta(1))*obj.pde.stiffMat;
      rhs = -obj.M0.stiffMat*(u0*reshape(obj.data.alpha(end:-1:2),[],1)) + ...
          (diff(I)*sum(obj.data.beta))*obj.pde.loadVec;
      if numel(obj.data.beta)>1
        rhs = rhs - diff(I)*(obj.pde.stiffMat*(u0*reshape(obj.data.beta(end:-1:2),[],1)));
      end
      R = obj.solver.solve(lhs, rhs, obj.freeI, obj.freeJ, obj.pde.getShift());
    end
  end
end