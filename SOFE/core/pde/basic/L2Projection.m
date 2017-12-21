classdef L2Projection < PDE
% Id(1) = FId(f)
  methods % constructor
    function obj = L2Projection(data, fesTrial, varargin)
      if ~isempty(varargin), fesTest = varargin{:}; else, fesTest = fesTrial; end
      opList = {OpIdId(1, 0, fesTrial, varargin{:}), FcId(data.f, fesTest, 0)};
      lhs.sys = {{1}}; rhs.sys = {{2}};
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end