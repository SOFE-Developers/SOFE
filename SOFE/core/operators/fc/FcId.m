classdef FcId < Functional % ( F, V )
  methods % constructor
    function obj = FcId(coeff, fes, codim, varargin) % [loc]
      obj = obj@Functional(coeff, fes, codim, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      basis = [];
      if obj.codim  <  obj.fes.element.dimension
        basis = obj.fes.evalGlobalBasis([], obj.codim, 0, {k}); % nExnBxnPxnC
      end
      R = integrate(obj, basis, k);
    end
  end
end
