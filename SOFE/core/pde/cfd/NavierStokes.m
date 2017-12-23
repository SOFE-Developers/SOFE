classdef NavierStokes < PDE
  methods % constructor
    function obj = NavierStokes(data, fesV, fesP)
      opList = {OpSGradSGrad(data.nu, fesV), ...
                OpGradId(@(x,t,U)U{1}, fesV, fesV), ...
                OpDivId(1, fesV, fesP), ...
                OpIntId(fesP)};
      lhs.sys = {{1,2} {3} {}; {3} {} {4}; {} {4} {}};
      lhs.coeff = {{} {-1.0} {}; {} {} {}; {} {} {}};
      lhs.adj = {{} {1} {}; {} {} {1}; {} {} {}};
      rhs.sys = {{}; {}; {}};
      %
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end
