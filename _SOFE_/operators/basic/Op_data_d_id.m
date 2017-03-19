classdef Op_data_d_id < Operator % ( c*dx_k(u), v )
  properties
    kIdx
  end
  methods % constructor
    function obj = Op_data_d_id(coeff, kIdx, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
      obj.kIdx = kIdx;
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnD
      dBasisJ = dBasisJ(:,:,:,obj.kIdx); % nExnBxnP
      basisI = obj.feSpaceTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnP
      R = obj.integrate(false, basisI, dBasisJ, k);
    end
  end
end
