classdef OpCurlCurl < Operator % ( c*Curl(U), Curl(V) )
  methods % constructor
    function obj = OpCurlCurl(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      if obj.fesTrial.element.dimension == 2
        dBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
        dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
        try
          curlBasisI = dBasisI(:,:,:,2,1) - dBasisI(:,:,:,1,2); % nExnBxnP
          curlBasisJ = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
        catch
          curlBasisI(:,:,:,1) = dBasisI(:,:,:,1,2); % nExnBxnP
          curlBasisI(:,:,:,2) = -dBasisI(:,:,:,1,1); % nExnBxnP
          curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,1,2); % nExnBxnP
          curlBasisJ(:,:,:,2) = -dBasisJ(:,:,:,1,1); % nExnBxnP
        end
        R = obj.integrate(curlBasisI, curlBasisJ, k);
      else
%         dBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
%         curlBasisI(:,:,:,1) = dBasisI(:,:,:,3,2) - dBasisI(:,:,:,2,3); % nExnBxnP
%         curlBasisI(:,:,:,2) = dBasisI(:,:,:,1,3) - dBasisI(:,:,:,3,1); % nExnBxnP
%         curlBasisI(:,:,:,3) = dBasisI(:,:,:,2,1) - dBasisI(:,:,:,1,2); % nExnBxnP
        curlBasisI = obj.fesTest.evalGlobalBasis([], 0, 'curl', {k}); % nExnBxnPxnC
        if obj.fesTrial ~= obj.fesTest
%           dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
%           curlBasisJ(:,:,:,1) = dBasisJ(:,:,:,3,2) - dBasisJ(:,:,:,2,3); % nExnBxnP
%           curlBasisJ(:,:,:,2) = dBasisJ(:,:,:,1,3) - dBasisJ(:,:,:,3,1); % nExnBxnP
%           curlBasisJ(:,:,:,3) = dBasisJ(:,:,:,2,1) - dBasisJ(:,:,:,1,2); % nExnBxnP
          curlBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 'curl', {k}); % nExnBxnPxnC
          R = obj.integrate(curlBasisI, curlBasisJ, k);
        else
          R = obj.integrate(curlBasisI, curlBasisI, k);
        end
      end
    end
    function R = getScaling(obj, nRef)
      switch obj.fesTrial.element.conformity
        case {'H1', 'L2'}
          R = 2^((nRef*(2-obj.fesTrial.element.dimension)));
        case 'HDiv'
          assert(false, 'not allowed');
        case 'HRot'
          R = 2^((nRef*(4-obj.fesTrial.element.dimension)));
      end
    end
  end
end
