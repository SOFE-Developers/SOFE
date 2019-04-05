classdef L2Projection < PDE
  methods % constructor
    function obj = L2Projection(data, fesTrial, varargin)
      if ~isempty(varargin), fesTest = varargin{:}; else, fesTest = fesTrial; end
      try coef = data.c; catch, coef = 1.0; end
      opList = {OpIdId(coef, 0, fesTrial, varargin{:}), FcId(data.f, fesTest, 0)};
%       opList = {OpMatIdId(coef, 0, fesTrial, varargin{:}), FcId(data.f, fesTest, 0)};
      lhs.sys = {{1}}; rhs.sys = {{2}};
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end