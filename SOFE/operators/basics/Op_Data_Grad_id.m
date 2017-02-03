classdef Op_Data_Grad_id < Op_Data_GRAD_Id % ( C*Grad(u), v )
  methods % constructor
    function obj = Op_Data_Grad_id(coeff, feSpaceTrial, varargin)
      obj = obj@Op_Data_GRAD_Id(coeff, feSpaceTrial, varargin{:});
    end
  end
end
