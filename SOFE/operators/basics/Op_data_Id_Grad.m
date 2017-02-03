classdef Op_data_Id_Grad < Operator % ( c*U, grad(v) )
  methods % constructor
    function obj = Op_data_Id_Grad(coeff, feSpaceTrial, feSpaceTest)
      obj = obj@Operator(coeff, feSpaceTrial, feSpaceTest);
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      basisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 0, varargin{:}); % nExnBxnPxnCx1
      gradBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPx1xnD
      R = obj.integrate(true, gradBasisI, basisJ, varargin{:});
    end
  end
end
