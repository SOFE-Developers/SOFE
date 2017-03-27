classdef Op_data_Curl_Curl < Operator % ( c*Curl(U), Curl(V) )
  methods % constructor
    function obj = Op_data_Curl_Curl(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k});
      dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      if obj.fesTrial.element.dimension == 2
        curlBasisI = dBasisI(:,:,:,2,1) - dBasisI(:,:,:,1,2); % nExnBxnP
        curlBasisJ = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
      else
        curlBasisI(:,:,:,1) = dBasisI(:,:,:,3,2) - dBasisI(:,:,:,2,3); % nExnBxnP
        curlBasisI(:,:,:,2) = dBasisI(:,:,:,1,3) - dBasisI(:,:,:,3,1); % nExnBxnP
        curlBasisI(:,:,:,3) = dBasisI(:,:,:,2,1) - dBasisI(:,:,:,1,2); % nExnBxnP
        curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,3,2) - dBasisJ(:,:,:,2,3); % nExnBxnP
        curlBasisJ(:,:,:,2) = dBasisJ(:,:,:,1,3) - dBasisJ(:,:,:,3,1); % nExnBxnP
        curlBasisJ(:,:,:,3) = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
      end
      R = obj.integrate(true, curlBasisI, curlBasisJ, k);
    end
  end
end
