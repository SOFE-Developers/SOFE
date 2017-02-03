classdef Op_data_Grad_Grad < Op_data_GRAD_GRAD % ( c*Grad(u), Grad(v) )
  methods % constructor
    function obj = Op_data_Grad_Grad(coeff, feSpaceTrial, varargin)
      obj = obj@Op_data_GRAD_GRAD(coeff, feSpaceTrial, varargin{:});
    end
  end
end
