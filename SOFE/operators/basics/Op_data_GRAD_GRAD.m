classdef Op_data_GRAD_GRAD < Operator % ( c*grad(U), grad(V) )
  methods % constructor
    function obj = Op_data_GRAD_GRAD(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      gradBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, k); % nExnBxnPxnCxnD
      if obj.isGalerkin
        gradBasisI = gradBasisJ;
      else
        gradBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, k);
      end
      R = obj.integrate(true, gradBasisI, gradBasisJ, k);
    end
  end
end
