classdef FcDiv < Functional % ( F, div(V) )
  methods % constructor
    function obj = FcDiv(coeff, fes, varargin) % [loc]
      obj = obj@Functional(coeff, fes, 0, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      basis = [];
      if obj.codim  <  obj.fes.element.dimension
        basis = obj.fes.evalGlobalBasis([], obj.codim, 1, {k}); % nExnBxnPxnCxnD
        basis = basis(:,:,:,1,1) + basis(:,:,:,2,2);
      end
      R = integrate(obj, basis, k);
    end
  end
end
