classdef Poisson2 < PDE2
% GradGrad(a) + Id_Gamma(h) = FId(f) + FId_Gamma(g)
  methods % constructor
    function obj = Poisson2(data, fesTrial, varargin)
      opList = {OpGradGrad(data.a, fesTrial, varargin{:}), ...
                OpGradId([1 1], fesTrial, varargin{:})};
      lhs = {{1,2}}; opFactor = {{1.0, 10.0}}; opAdj = {{0 0}};
      %
      fcList = {FcId(data.f, opList{1}.fesTest, 0)};
      rhs = {{1}}; fcFactor = {{1}};
      %
      obj = obj@PDE2(lhs, rhs, opList, opFactor, opAdj, fcList, fcFactor);
    end
  end
end