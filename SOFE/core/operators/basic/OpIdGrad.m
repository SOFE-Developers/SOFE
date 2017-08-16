classdef OpIdGrad < Operator % (c*V, Grad(u))
  methods % constructor
    function obj = OpIdGrad(coeff, fesTrial, fesTest)
      obj = obj@Operator(coeff, fesTrial, fesTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      basisJ = obj.fesTrial.evalGlobalBasis([], 0, 0, {k}); % nExnBxnPxnCx1
      gradBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPx1xnD
      R = obj.integrate(true, gradBasisI, basisJ, k);
    end
  end
end
