classdef Theta < TimeStep
  properties
    freeI
    freeJ
    theta
  end
  methods % constructor
    function obj = Theta(theta, M0, pde, solver)
      obj = obj@TimeStep([], M0, pde, solver);
      obj.nS = 1;
      obj.theta = theta;
      obj.setFreeDoFs();
    end
    function setFreeDoFs(obj)
      [obj.freeI, obj.freeJ] = obj.pde.getFreeDoFs();
    end
  end
  methods % integrate
    function R = compute(obj, I, u0)
      obj.pde.setState(I(1), u0); obj.pde.assemble();
      rhs = obj.M0.stiffMat*u0 + ((1-obj.theta)*diff(I))*(obj.pde.loadVec - obj.pde.stiffMat*u0);
      obj.pde.setState(I(2), u0); obj.pde.assemble();
      lhs = obj.M0.stiffMat + (obj.theta*diff(I))*obj.pde.stiffMat;
      rhs = rhs + (obj.theta*diff(I))*obj.pde.loadVec;
      R = obj.solver.solve(lhs, rhs, obj.freeI, obj.freeJ, obj.pde.getShift());
    end
  end
end