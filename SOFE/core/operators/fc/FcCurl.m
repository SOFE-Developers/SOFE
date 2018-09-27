classdef FcCurl < Functional % ( F, curl(V) )
  methods % constructor
    function obj = FcCurl(coeff, fes, varargin) % [loc]
      obj = obj@Functional(coeff, fes, 0, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      curlBasis = [];
      if obj.codim  <  obj.fes.element.dimension
        basis = obj.fes.evalGlobalBasis([], obj.codim, 1, {k}); % nExnBxnPxnCxnD
        if obj.fes.element.dimension == 2
          curlBasis = basis(:,:,:,2,1) - basis(:,:,:,1,2);
        else
          curlBasis(:,:,:,3) = basis(:,:,:,2,1) - basis(:,:,:,1,2); % nExnBxnP
          curlBasis(:,:,:,2) = basis(:,:,:,1,3) - basis(:,:,:,3,1); % nExnBxnP
          curlBasis(:,:,:,1) = basis(:,:,:,3,2) - basis(:,:,:,2,3); % nExnBxnP
        end
      end
      R = integrate(obj, curlBasis, k);
    end
  end
end
