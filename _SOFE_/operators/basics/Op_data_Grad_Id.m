classdef Op_data_Grad_Id < Operator % (c*Grad(u), V )
  methods % constructor
    function obj = Op_data_Grad_Id(coeff, feSpaceTrial, feSpaceTest)
      obj = obj@Operator(coeff, feSpaceTrial, feSpaceTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      gradBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, k); % nExnBxnPx1xnD
      basisI = obj.feSpaceTest.evalGlobalBasis([], 0, 0, k); % nExnBxnPxnCx1
      R = obj.integrate(true, basisI, gradBasisJ, k);
    end
  end
end
