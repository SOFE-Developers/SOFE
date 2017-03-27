classdef Op_DATA_Grad_Grad < Op_DATA_GRAD_GRAD % ( CC*Grad(u), Grad(v) )
  methods % constructor
    function obj = Op_DATA_Grad_Grad(coeff, fesTrial, varargin)
      obj = obj@Op_DATA_GRAD_GRAD(coeff, fesTrial, varargin{:});
    end
  end
end
