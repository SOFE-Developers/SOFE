classdef Op_data_SYMGRAD_SYMGRAD < Operator % ( c*SYMGRAD(U), SYMGRAD(V) )
  methods % constructor
    function obj = Op_data_SYMGRAD_SYMGRAD(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      sGradBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      sGradBasisJ = 0.5*(sGradBasisJ + permute(sGradBasisJ, [1 2 3 5 4]));
      if obj.isGalerkin
        sGradBasisI = sGradBasisJ;
      else
        sGradBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, {k});
        sGradBasisI = 0.5*(sGradBasisI + permute(sGradBasisI, [1 2 3 5 4]));
      end
      R = obj.integrate(true, sGradBasisI, sGradBasisJ, k);
    end
  end
end
