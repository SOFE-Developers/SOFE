classdef Op_data_GRAD_GRAD < Operator % ( c*GRAD(U), GRAD(V) )
  methods % constructor
    function obj = Op_data_GRAD_GRAD(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      gradBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      gradBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k});
      R = obj.integrate(true, gradBasisI, gradBasisJ, k);
    end
  end
end
