classdef OpCurlId < Operator % ( c*Curl(U), V )
  methods % constructor
    function obj = OpCurlId(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      if obj.fesTrial.element.dimension == 2
        dBasis = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
        try
          curlBasis = dBasis(:,:,:,2,1) - dBasis(:,:,:,1,2); % nExnBxnP
        catch
          curlBasis(:,:,:,1) = dBasis(:,:,:,1,2); % nExnBxnP
          curlBasis(:,:,:,2) = -dBasis(:,:,:,1,1); % nExnBxnP
        end
      else
%         dBasis = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
%         curlBasis(:,:,:,1) = dBasis(:,:,:,3,2) - dBasis(:,:,:,2,3); % nExnBxnP
%         curlBasis(:,:,:,2) = dBasis(:,:,:,1,3) - dBasis(:,:,:,3,1); % nExnBxnP
%         curlBasis(:,:,:,3) = dBasis(:,:,:,2,1) - dBasis(:,:,:,1,2); % nExnBxnP
        curlBasis = obj.fesTrial.evalGlobalBasis([], 0, 'curl', {k}); % nExnBxnPxnC
      end
      %
      basis = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnPxnC
      R = obj.integrate(basis, curlBasis, k);
    end
    function R = getScaling(obj, nRef)
      R = 2^(nRef*0*obj.fesTrial.element.dimension);
    end
  end
end
