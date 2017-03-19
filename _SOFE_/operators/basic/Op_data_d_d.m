classdef Op_data_d_d < Operator % ( c*dx_k(u), dx_l(v) )
  properties
    kIdx
    lIdx
  end
  methods % constructor
    function obj = Op_data_d_d(coeff, kIdx, lIdx, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
      obj.kIdx = kIdx;
      obj.lIdx = lIdx;
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnD
      dBasisJ = dBasisJ(:,:,:,obj.kIdx); % nExnBxnP
      dBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnD
      dBasisI = dBasisI(:,:,:,obj.lIdx); % nExnBxnP
      R = obj.integrate(false, basisI, dBasisJ, k);
    end
  end
end
