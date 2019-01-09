classdef Stokes < PDE
  methods % constructor
    function obj = Stokes(data, fesV, fesP)
      opList = {OpGradGrad(data.nu, fesV), ...
                OpDivId(-1.0, fesV, fesP), ...
                FcId(data.f, fesV, 0)};
      lhs.sys = {{1} {2}; ...
                 {2} {}};
      lhs.adj = {{} {1}; ...
                 {} {}};
      rhs.sys = {{3}; {}};
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end
