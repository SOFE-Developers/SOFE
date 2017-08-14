classdef OpGradId < Operator % (c*Grad(u), V )
  methods % constructor
    function obj = OpGradId(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      basisI = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnPxnCx1
      dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPx1xnD
      hasCoeff = true;
      if obj.fesTest.element.getNC() == obj.fesTrial.element.getNC()
        if isnumeric(obj.data)
          C = obj.fesTrial.evalDoFVector(obj.data, [], 0, 0, {k}); % nExnPx(nD*nD)
        else
          C = obj.fesTrial.evalFunction(obj.data, [], 0, obj.state, obj.dState, {k}); % nExnPxnW
        end
        C = permute(C, [1 4 2 5 3]); % nEx1xnPx1xnW
        dBasisJ = sum(bsxfun(@times, dBasisJ, C), 5); % nExnBxnPxnC
        hasCoeff = false;
      end
      R = obj.integrate(hasCoeff, basisI, dBasisJ, k);
    end
  end
end
