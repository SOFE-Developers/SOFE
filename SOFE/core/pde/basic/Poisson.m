classdef Poisson < PDE
  methods % constructor
    function obj = Poisson(data, fesTrial, varargin)
      if ~isempty(varargin), fesTest = varargin{:}; else, fesTest = fesTrial; end
      opList = {OpGradGrad(data.a, fesTrial, varargin{:}), ...
                FcId(data.f, fesTest, 0)};
      lhs.sys = {{1}};
      rhs.sys = {{2}};
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end