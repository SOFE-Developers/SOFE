classdef LotkaVolterra < PDE
  methods % constructor
    function obj = LotkaVolterra(data, fes)
      lap = Op_data_GRAD_GRAD(data.a{1}, fes, fes);
      f1 = Fc_Data_Id(data.f{1}, fes, 0);
      f2 = Fc_Data_Id(data.f{2}, fes, 0);
      %
      lhs = {{lap},  {}; ...
              {}, {lap}};
      rhs = {{f1}; {f2}};
      obj = obj@PDE(lhs, rhs);
    end
  end
end
