classdef Op_data_Grad_Id < Operator % (c*Grad(u), V )
  methods % constructor
    function obj = Op_data_Grad_Id(coeff, fesTrial, fesTest)
      obj = obj@Operator(coeff, fesTrial, fesTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      gradBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPx1xnD
      basisI = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnPxnCx1
      R = obj.integrate(true, basisI, gradBasisJ, k);
    end
  end
end
