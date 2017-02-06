classdef Op_data_id_div < Operator % (c*u, div(V) )
  methods % constructor
    function obj = Op_data_id_div(coeff, feSpaceTrial, feSpaceTest)
      obj = obj@Operator(coeff, feSpaceTrial, feSpaceTest);
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      basisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 0, varargin{:}); % nExnBxnP
      dBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnCxnD
      nD = size(dBasisI,5);
      divBasisI = dBasisI(:,:,:,1,1); % nExnBxnP
      for d = 2:nD
        divBasisI = divBasisI + dBasisI(:,:,:,d,d); % nExnBxnP
      end
      R = obj.integrate(true, divBasisI, basisJ, varargin{:});
    end
  end
end
