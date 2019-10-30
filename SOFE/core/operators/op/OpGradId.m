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
      if obj.fesTest.element.getNC() == obj.fesTrial.element.getNC()
        if isnumeric(obj.data)
          if size(obj.data,1)==1
            C = permute(obj.data, [3 1 2]); % 1x1xnC
          else
            C = obj.fesTrial.evalDoFVector(obj.data, [], 0, 0, {k}); % nExnPxnC
          end
        else
          S = obj.observers{1}.evalState(k);
          C = obj.fesTrial.evalFunction(obj.data, [], 0, S, {k}); % nExnPxnC
        end
        C = permute(C, [1 4 2 5 3]); % nEx1xnPx1xnW
        dBasisJ = sum(bsxfun(@times, dBasisJ, C), 5); % nExnBxnPxnC
        obj.hasCoeff = false;
      end
      R = obj.integrate(basisI, dBasisJ, k);
    end
    function R = getScaling(obj, nRef)
%       warning('TODO: scaling!');
%       if obj.fesTest.element.getNC() ~= obj.fesTrial.element.getNC()
%         R = 2^((nRef*0*obj.fesTrial.element.dimension));
%       else
        R = 2^((nRef*(1-obj.fesTrial.element.dimension)));
%       end
    end
  end
end
