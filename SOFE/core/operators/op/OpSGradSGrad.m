classdef OpSGradSGrad < Operator % ( c*SGRAD(U), SGRAD(V) )
  methods % constructor
    function obj = OpSGradSGrad(coeff, fes)
      obj = obj@Operator(coeff, fes);
    end
  end
  methods
    function R = assembleOp(obj, k)
      sGradBasis = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      sGradBasis = 0.5*(sGradBasis + permute(sGradBasis, [1 2 3 5 4]));
      R = obj.integrate(sGradBasis, sGradBasis, k);
    end
    function R = getScaling(obj, nRef)
      R = 2^((nRef*(2-obj.fesTrial.element.dimension)));
    end
  end
end
