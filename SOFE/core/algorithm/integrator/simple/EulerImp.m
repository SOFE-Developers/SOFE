classdef EulerImp < TimeStep
  methods % constructor
    function obj = EulerImp(M0, pde, solver)
      obj = obj@TimeStep([], M0, pde, solver);
      obj.nS = 1; obj.nK = 1;
    end
  end
  methods % integrate
    function R = compute(obj, I, u0)
      obj.pde.setState(I(2), u0);
      obj.pde.assemble();
      lhs = obj.M0.stiffMat + diff(I)*obj.pde.stiffMat;
      rhs = diff(I)*obj.pde.loadVec + obj.M0.stiffMat*u0;
      R = obj.solver.solve(lhs, rhs, obj.freeI, obj.freeJ, obj.pde.getShift());
    end
  end
end