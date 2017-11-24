classdef OpCurlId < Operator % ( c*Curl(U), V )
  methods % constructor
    function obj = OpCurlId(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      if obj.fesTrial.element.dimension == 2
        try
          curlBasisJ = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
        catch
          curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,1,2); % nExnBxnP
          curlBasisJ(:,:,:,2) = -dBasisJ(:,:,:,1,1); % nExnBxnP
        end
      else
        curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,3,2) - dBasisJ(:,:,:,2,3); % nExnBxnP
        curlBasisJ(:,:,:,2) = dBasisJ(:,:,:,1,3) - dBasisJ(:,:,:,3,1); % nExnBxnP
        curlBasisJ(:,:,:,3) = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
      end
      %
      basis = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnPxnC
      R = obj.integrate(true, basis, curlBasisJ, k);
    end
  end
end
