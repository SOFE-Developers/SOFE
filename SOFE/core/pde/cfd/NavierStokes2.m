classdef NavierStokes2 < PDE
  methods % constructor
    function obj = NavierStokes2(data, fesV, fesP)
      opList = {OpGradGrad(data.nu, fesV), ...
                OpGradId(@(x,t,U)U{1}, fesV, fesV), ...
                OpDivId(-1.0, fesV, fesP), ...
                FcId(data.f, fesV, 0)};
      try
        opList{5} = OpGradGrad(data.h, fesP);
      catch
        opList{5} = OpGradGrad(0, fesP);
      end
      lhs.sys = {{1,2} {3}; ...
                   {3} {5}};
      lhs.coeff = {{} {}; ...
                   {} {}};
      lhs.adj = {{} {1}; ...
                 {} {}};
      rhs.sys = {{4}; {}};
      %
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end
