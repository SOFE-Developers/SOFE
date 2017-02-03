classdef Op_data_Curl_Curl < Operator % ( c*curl(U), curl(V) )
  methods % constructor
    function obj = Op_data_Curl_Curl(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnCxnD
      if obj.feSpaceTrial.element.dimension == 2
        curlBasisJ = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
      else
        curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,3,2) - dBasisJ(:,:,:,2,3); % nExnBxnP
        curlBasisJ(:,:,:,2) = dBasisJ(:,:,:,1,3) - dBasisJ(:,:,:,3,1); % nExnBxnP
        curlBasisJ(:,:,:,3) = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
      end
      if obj.isGalerkin
        curlBasisI = curlBasisJ;
      else
        dBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, varargin{:});
        if obj.feSpaceTrial.element.dimension == 2
          curlBasisI = dBasisI(:,:,:,2,1) - dBasisI(:,:,:,1,2); % nExnBxnP
        else
          curlBasisI(:,:,:,1) = dBasisI(:,:,:,3,2) - dBasisI(:,:,:,2,3); % nExnBxnP
          curlBasisI(:,:,:,2) = dBasisI(:,:,:,1,3) - dBasisI(:,:,:,3,1); % nExnBxnP
          curlBasisI(:,:,:,3) = dBasisI(:,:,:,2,1) - dBasisI(:,:,:,1,2); % nExnBxnP
        end
      end
      R = obj.integrate(true, curlBasisI, curlBasisJ, varargin{:});
    end
  end
end
