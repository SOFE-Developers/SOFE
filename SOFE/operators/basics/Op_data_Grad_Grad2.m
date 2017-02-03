classdef Op_data_Grad_Grad2 < Op_data_GRAD_GRAD % ( c*grad(u), grad(v) )
  methods % constructor
    function obj = Op_data_Grad_Grad2(coeff, feSpaceTrial, varargin)
      obj = obj@Op_data_GRAD_GRAD(coeff, feSpaceTrial, varargin{:});
    end
  end
end
