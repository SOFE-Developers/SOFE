classdef Op_Data_GRAD_Id2 < Operator % ( C*grad(u), v )
  methods % constructor
    function obj = Op_Data_GRAD_Id2(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      C = obj.feSpaceTrial.evalFunction(obj.data, [], 0, obj.state, varargin{:}); % nExnPxnW
      C = permute(C, [1 4 2 5 3]); % nEx1xnPx1xnW
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnCxnW
      dBasisJ = sum(bsxfun(@times, dBasisJ, C), 5); % nExnBxnPxnC
      basisI = obj.feSpaceTest.evalGlobalBasis([], 0, 0, varargin{:}); % nExnBxnPxnC
      R = obj.integrate(false, basisI, dBasisJ, varargin{:});
    end
  end
end
