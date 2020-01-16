classdef OpGradGrad < Operator % ( c*GRAD(U), GRAD(V) )
  methods % constructor
    function obj = OpGradGrad(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      gradBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      if obj.fesTrial ~= obj.fesTest
        gradBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k});
        R = obj.integrate(gradBasisI, gradBasisJ, k);
      else
        R = obj.integrate(gradBasisJ, gradBasisJ, k);
      end
    end
    function R = getScaling(obj, nRef)
      R = 2^((nRef*(2-obj.fesTrial.element.dimension)));
    end
  end
end
