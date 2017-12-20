classdef Poisson2 < PDE2
  methods % constructor
    function obj = Poisson2(data, fesTrial, varargin)
      if ~isempty(varargin), fesTest = varargin{:}; else, fesTest = fesTrial; end
      opList = {OpGradGrad(data.a, fesTrial, fesTest), ...
                FcId(data.f, fesTest, 0)};
      lhs.sys = {{1}};
      rhs.sys = {{2}};
      obj = obj@PDE2(opList, lhs, rhs);
    end
  end
end