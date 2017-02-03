classdef Op_data_d_d < Operator % ( c*d_dx_k(u), d_dx_k(v) )
  properties
    kIndex
    lIndex
  end
  methods % constructor
    function obj = Op_data_d_d(coeff, kIndex, lIndex, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
      obj.kIndex = kIndex;
      obj.lIndex = lIndex;
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnD
      dBasisJ = dBasisJ(:,:,:,obj.kIndex); % nExnBxnP
      dBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnD
      dBasisI = dBasisI(:,:,:,obj.lIndex); % nExnBxnP
      R = obj.integrate(false, basisI, dBasisJ, varargin{:});
    end
  end
end
