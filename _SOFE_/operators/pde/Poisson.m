classdef Poisson < PDE
% GradGrad(a) + Id_Gamma(h) = FId(f) + FId_Gamma(g)
  methods % constructor
    function obj = Poisson(data, fesTrial, varargin)
      lhs = {Op_data_GRAD_GRAD(data.a, fesTrial, varargin{:})};
      rhs = {Fc_Data_Id(data.f, lhs{1}.fesTest, 0)};
      try % Neumann BC
        rhs{2} = Fc_Data_Id(data.g, lhs{1}.fesTest, 1);
      end
      try % Robin BC
        lhs{2} = Op_data_Id_Id(data.h, 1, lhs{1}.fesTrial, lhs{1}.fesTest);
      end
      obj = obj@PDE({lhs}, {rhs});
    end
  end
end