classdef L2Projection < PDE
% Id(1) = FId(f)
  methods % constructor
    function obj = L2Projection(data, fesTrial, varargin)
      lhs = {Op_data_Id_Id(1, 0, fesTrial, varargin{:})};
      rhs = {Fc_Data_Id(data.f, lhs{1}.feSpaceTest, 0)};
      obj = obj@PDE({lhs}, {rhs});
    end
  end
end