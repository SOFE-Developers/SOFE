classdef CDR < PDE
% GradGrad(a) + Grad(b) + Id(c) + Id_Gamma(h) = FId(f) + FId_Gamma(g)
  methods % constructor
    function obj = CDR(data, fesTrial, varargin)
      lap = Op_data_GRAD_GRAD(data.a, fesTrial, varargin{:});
      lhs = {lap};
      try
        lhs = [lhs, {Op_Data_Grad_id(data.b, fesTrial, varargin{:})}];
      end
      try
        lhs = [lhs, {Op_data_Id_Id(data.c, 0, fesTrial, varargin{:})}];
      end
      try
        lhs = [lhs, {Op_data_Id_Id(data.h, 1, lap.fesTrial, lap.fesTest)}];
      end
      %
      try
        rhs = {Fc_Data_Id(data.f, lap.fesTest, 0)};
      catch
        rhs = {Fc_Data_Id(@(x)0*x(:,1), lap.fesTest, 0)};
      end
      try
        rhs = [rhs,{Fc_Data_Id(data.g, lap.fesTest, 1)}];
      end
      obj = obj@PDE({lhs}, {rhs});
    end
  end
end
