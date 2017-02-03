classdef Op_data_Curl_Id < Operator % ( c*curl(U), V )
  methods % constructor
    function obj = Op_data_Curl_Id(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, varargin)
      dBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, varargin{:}); % nExnBxnPxnCxnD
      curlBasisJ = zeros([size(dBasisJ,1),size(dBasisJ,2),size(dBasisJ,3),size(dBasisJ,4)]); % nExnBxnPxnC
      curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,3,2)-dBasisJ(:,:,:,2,3);
      curlBasisJ(:,:,:,2) = dBasisJ(:,:,:,1,3)-dBasisJ(:,:,:,3,1);
      curlBasisJ(:,:,:,3) = dBasisJ(:,:,:,2,1)-dBasisJ(:,:,:,1,2);
      basis = obj.feSpaceTest.evalGlobalBasis([], 0, 0, varargin{:}); % nExnBxnPxnC
      R = obj.integrate(true, basis, curlBasisJ, varargin{:});
    end
  end
end
