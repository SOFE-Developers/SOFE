classdef Op_data_d_id < Operator % ( c*d_dx_k(u), v )
  properties
    kIndex
  end
  methods % constructor
    function obj = Op_data_d_id(coeff, kIndex, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
      obj.kIndex = kIndex;
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnD
      dBasisJ = dBasisJ(:,:,:,obj.kIndex); % nExnBxnP
      basisI = obj.feSpaceTest.evalGlobalBasis([], 0, 0, varargin{:}); % nExnBxnP
      R = obj.integrate(false, basisI, dBasisJ, varargin{:});
    end
  end
end
