classdef OpDId < Operator % ( c*dx_k(u), v )
  properties
    kIdx
  end
  methods % constructor
    function obj = OpDId(coeff, kIdx, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
      obj.kIdx = kIdx;
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnD
      dBasisJ = dBasisJ(:,:,:,:,obj.kIdx); % nExnBxnPxnCxnD
      basisI = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnP
      R = obj.integrate(basisI, dBasisJ, k);
    end
  end
end
