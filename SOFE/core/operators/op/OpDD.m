classdef OpDD < Operator % ( c*dx_k(u), dx_l(v) )
  properties
    kIdx
    lIdx
  end
  methods % constructor
    function obj = OpDD(coeff, kIdx, lIdx, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
      obj.kIdx = kIdx;
      obj.lIdx = lIdx;
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnD
      dBasisJ = dBasisJ(:,:,:,obj.kIdx); % nExnBxnP
      dBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnD
      dBasisI = dBasisI(:,:,:,obj.lIdx); % nExnBxnP
      R = obj.integrate(dBasisI, dBasisJ, k);
    end
  end
end
