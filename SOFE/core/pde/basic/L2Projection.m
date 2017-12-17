classdef L2Projection < PDE
% Id(1) = FId(f)
  methods % constructor
    function obj = L2Projection(data, fesTrial, varargin)
      lhs = {OpIdId(1, 0, fesTrial, varargin{:})};
      rhs = {FcId(data.f, lhs{1}.fesTest, 0)};
      obj = obj@PDE({lhs}, {rhs});
    end
  end
end