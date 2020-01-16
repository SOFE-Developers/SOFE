classdef FcDiv < Functional % ( F, div(V) )
  methods % constructor
    function obj = FcDiv(coeff, fes, varargin) % [loc]
      obj = obj@Functional(coeff, fes, 0, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      divBasis = [];
      if obj.codim  <  obj.fes.element.dimension
%       dBasis = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
%       nC = size(dBasis,4);
%       divBasis = dBasis(:,:,:,1,1); % nExnBxnP
%       for d = 2:nC
%         divBasis = divBasis + dBasis(:,:,:,d,d); % nExnBxnP
%       end
        divBasis = obj.fes.evalGlobalBasis([], obj.codim, 'div', {k}); % nExnBxnPxnCxnD
      end
      R = integrate(obj, divBasis, k);
    end
  end
end
