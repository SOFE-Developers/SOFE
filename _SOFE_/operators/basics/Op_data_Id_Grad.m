classdef Op_data_Id_Grad < Operator % (c*Grad(u), V )
  methods % constructor
    function obj = Op_data_Id_Grad(coeff, feSpaceTrial, feSpaceTest)
      obj = obj@Operator(coeff, feSpaceTrial, feSpaceTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      basisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 0, k); % nExnBxnPxnCx1
      gradBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, k); % nExnBxnPx1xnD
      R = obj.integrate(true, gradBasisI, basisJ, k);
    end
  end
end
