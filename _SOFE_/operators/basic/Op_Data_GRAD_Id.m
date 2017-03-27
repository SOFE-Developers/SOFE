classdef Op_Data_GRAD_Id < Operator % ( C*GRAD(U), V )
  methods % constructor
    function obj = Op_Data_GRAD_Id(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      C = obj.fesTrial.evalFunction(obj.data, [], 0, obj.state, {k}); % nExnPxnW
      C = permute(C, [1 4 2 5 3]); % nEx1xnPx1xnW
      dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnW
      dBasisJ = sum(bsxfun(@times, dBasisJ, C), 5); % nExnBxnPxnC
      basisI = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnPxnC
      R = obj.integrate(false, basisI, dBasisJ, k);
    end
  end
end
