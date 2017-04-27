classdef Op_Data_Grad_id < Op_Data_GRAD_Id % ( C*Grad(u), v )
  methods % constructor
    function obj = Op_Data_Grad_id(coeff, fesTrial, varargin)
      obj = obj@Op_Data_GRAD_Id(coeff, fesTrial, varargin{:});
    end
  end
end
