classdef FcDk < Functional % ( F, dk(v) )
  properties
    k
  end
  methods % constructor
    function obj = FcDk(coeff, k, fes, varargin) % [loc]
      obj = obj@Functional(coeff, fes, 0, varargin{:});
      obj.k = k;
    end
  end
  methods
    function R = assembleOp(obj, k)
      basis = [];
      if obj.codim  <  obj.fes.element.dimension
        basis = obj.fes.evalGlobalBasis([], obj.codim, 1, {k}); % nExnBxnPxnCxnD
      end
      R = integrate(obj, basis(:,:,:,:,obj.k), k);
    end
  end
end
