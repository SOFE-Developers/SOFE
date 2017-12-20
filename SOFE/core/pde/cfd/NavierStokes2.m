classdef NavierStokes2 < PDE2
  methods % constructor
    function obj = NavierStokes2(data, fesV, fesP)
      opList = {OpSGradSGrad(data.nu, fesV), ...
                OpGradId(@(x,t,U)U{1}, fesV, fesV), ...
                OpDivId(1, fesV, fesP)};
      lhs.sys = {{1,2} {3}; {3} {}};
      lhs.coeff = {{} {-1.0}; {} {}};
      lhs.adj = {{} {1}; {} {}};
      rhs.sys = {{}; {}};
      %
      obj = obj@PDE2(opList, lhs, rhs);
    end
  end
end
