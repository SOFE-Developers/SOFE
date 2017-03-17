classdef Op_data_Curl_Id < Operator % ( c*Curl(U), V )
  methods % constructor
    function obj = Op_data_Curl_Id(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      curlBasisJ = zeros(size(dBasisJ(:,:,:,:,1))); % nExnBxnPxnC
      curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,3,2)-dBasisJ(:,:,:,2,3);
      curlBasisJ(:,:,:,2) = dBasisJ(:,:,:,1,3)-dBasisJ(:,:,:,3,1);
      curlBasisJ(:,:,:,3) = dBasisJ(:,:,:,2,1)-dBasisJ(:,:,:,1,2);
      basis = obj.feSpaceTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnPxnC
      R = obj.integrate(true, basis, curlBasisJ, k);
    end
  end
end
