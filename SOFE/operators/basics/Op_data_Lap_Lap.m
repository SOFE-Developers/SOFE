classdef Op_data_Lap_Lap < Operator % ( a*laplace(U), laplace(V) )
  methods % constructor
    function obj = Op_data_Lap_Lap(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 2, varargin{:}); % nExnBxnPxnCxnDxnD
      [nE, nBJ, nP, nC, nD, ~] = size(dBasisJ);
      dBasisJ = reshape(dBasisJ, [nE,nBJ,nP,nC,nD*nD]);
      dBasisJ = sum(dBasisJ(:,:,:,:,1:(nD+1):(nD*nD)), 5); % nExnBxnPxnC
      if obj.isGalerkin
        dBasisI = dBasisJ;
      else
        dBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 2, varargin{:}); % nExnBxnPxnCxnDxnD
        nBI = size(dBasisI,2);
        dBasisI = reshape(dBasisI, [nE,nBI,nP,nC,nD*nD]);
        dBasisI = sum(dBasisI(:,:,:,:,1:(nD+1):(nD*nD)), 5); % nExnBxnPxnC
      end
      R = obj.integrate(true, dBasisI, dBasisJ, varargin{:});
    end
  end
end
