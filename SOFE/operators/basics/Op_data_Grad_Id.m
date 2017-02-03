classdef Op_data_Grad_Id < Operator % (c*grad(u), V )
  methods % constructor
    function obj = Op_data_Grad_Id(coeff, feSpaceTrial, feSpaceTest)
      obj = obj@Operator(coeff, feSpaceTrial, feSpaceTest);
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      gradBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPx1xnD
      basisI = obj.feSpaceTest.evalGlobalBasis([], 0, 0, varargin{:}); % nExnBxnPxnCx1
      R = obj.integrate(true, basisI, gradBasisJ, varargin{:});
    end
  end
end
