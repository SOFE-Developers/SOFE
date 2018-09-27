classdef FcGrad < Functional % ( F, grad(v) )
  methods % constructor
    function obj = FcGrad(coeff, fes, varargin) % [loc]
      obj = obj@Functional(coeff, fes, 0, varargin{:});
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
