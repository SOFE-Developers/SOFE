classdef CDR < PDE
% GradGrad(a) + Grad(b) + Id(c) + Id_Gamma(h) = FId(f) + FId_Gamma(g)
  methods % constructor
    function obj = CDR(data, fesTrial, varargin)
      lap = OpGradGrad(data.a, fesTrial, varargin{:});
      lhs = {lap};
      try
        lhs = [lhs, {OpGradId(data.b, fesTrial, varargin{:})}];
      end
      try
        lhs = [lhs, {OpIdId(data.c, 0, fesTrial, varargin{:})}];
      end
      try
        lhs = [lhs, {OpIdId(data.h, 1, lap.fesTrial, lap.fesTest)}];
      end
      %
      try
        rhs = {FcId(data.f, lap.fesTest, 0)};
      catch
        rhs = {FcId(@(x)0*x(:,1), lap.fesTest, 0)};
      end
      try
        rhs = [rhs,{FcId(data.g, lap.fesTest, 1)}];
      end
      obj = obj@PDE({lhs}, {rhs});
    end
  end
end
