classdef Op_data_SYMGRAD_SYMGRAD < Operator % ( c*symgrad(U), symgrad(V) )
  methods % constructor
    function obj = Op_data_SYMGRAD_SYMGRAD(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      sGradBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnCxnD
      sGradBasisJ = 0.5*(sGradBasisJ + permute(sGradBasisJ, [1 2 3 5 4]));
      if obj.isGalerkin
        sGradBasisI = sGradBasisJ;
      else
        sGradBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, varargin{:});
        sGradBasisI = 0.5*(sGradBasisI + permute(sGradBasisI, [1 2 3 5 4]));
      end
      R = obj.integrate(true, sGradBasisI, sGradBasisJ, varargin{:});
    end
  end
end
