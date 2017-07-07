classdef LotkaVolterra < PDE
  methods % constructor
    function obj = LotkaVolterra(data, fes)
      lap = OpGradGrad(data.a{1}, fes, fes);
      f1 = FcId(data.f{1}, fes, 0);
      f2 = FcId(data.f{2}, fes, 0);
      %
      lhs = {{lap},  {}; ...
              {}, {lap}};
      rhs = {{f1}; {f2}};
      obj = obj@PDE(lhs, rhs);
    end
  end
end
