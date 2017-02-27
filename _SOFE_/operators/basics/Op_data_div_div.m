classdef Op_data_div_div < Operator % ( c*div(U), div(V) )
  methods % constructor
    function obj = Op_data_div_div(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, k); % nExnBxnPxnCxnD
      nC = size(dBasisJ,4);
      divBasisJ = dBasisJ(:,:,:,1,1); % nExnBxnP
      for d = 2:nC
        divBasisJ = divBasisJ + dBasisJ(:,:,:,d,d); % nExnBxnP
      end
      if obj.isGalerkin
        divBasisI = divBasisJ;
      else
        dBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, k); % nExnBxnPxnCxnD
        divBasisI = dBasisI(:,:,:,1,1); % nExnBxnP
        for d = 2:nC
          divBasisI = divBasisI + dBasisI(:,:,:,d,d); % nExnBxnP
        end
      end
      R = obj.integrate(true, divBasisI, divBasisJ, k);
    end
  end
end
