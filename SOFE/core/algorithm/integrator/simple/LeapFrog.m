classdef LeapFrog < TimeStep
  properties
      freeI
      freeJ
  end
  methods % constructor
    function obj = LeapFrog(M0, pde, solver)
      obj = obj@TimeStep([], M0, pde, solver);
      obj.nS = 2;
      obj.setFreeDoFs();
    end
    function setFreeDoFs(obj)
      [obj.freeI, obj.freeJ] = obj.pde.getFreeDoFs();
    end
  end
  methods % integrate
    function R = compute(obj, I, u0)
      uneg1 = u0(:,1); u0 = u0(:,2);
      obj.pde.setState(I(2), u0);
      obj.pde.assemble();
      rhs = diff(I)^2*(obj.pde.loadVec - obj.pde.stiffMat*u0) + ...
            obj.M0.stiffMat*(2*u0 - uneg1);
      R = obj.solver.solve(obj.M0.stiffMat, rhs, obj.freeI, obj.freeJ, obj.pde.getShift());
    end
  end
end