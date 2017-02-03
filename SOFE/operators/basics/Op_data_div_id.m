classdef Op_data_div_id < Operator % (c*div(U), v )
  methods % constructor
    function obj = Op_data_div_id(coeff, feSpaceTrial, feSpaceTest)
      obj = obj@Operator(coeff, feSpaceTrial, feSpaceTest);
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnCxnD
      nD = size(dBasisJ,5);
      basisI = obj.feSpaceTest.evalGlobalBasis([], 0, 0, varargin{:}); % nExnBxnP
      divBasisJ = dBasisJ(:,:,:,1,1); % nExnBxnP
      for d = 2:nD
        divBasisJ = divBasisJ + dBasisJ(:,:,:,d,d); % nExnBxnP
      end
      R = obj.integrate(true, basisI, divBasisJ, varargin{:});
    end
  end
end
