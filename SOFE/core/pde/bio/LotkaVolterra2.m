classdef LotkaVolterra2 < PDE2
  methods % constructor
    function obj = LotkaVolterra2(data, fes)
      opList = {OpGradGrad(data.a{1}, fes, fes), ...
                FcId(data.f{1}, fes, 0), ...
                FcId(data.f{2}, fes, 0)};
      %
      lhs.sys = {{1}, {}; {} {1}};
      rhs.sys = {{2}; {3}};
      obj = obj@PDE2(opList, lhs, rhs);
    end
  end
end
