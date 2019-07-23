classdef OpDivDiv < Operator % ( c*div(U), div(V) )
  methods % constructor
    function obj = OpDivDiv(coeff, fes)
      obj = obj@Operator(coeff, fes);
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasis = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      nC = size(dBasis,4);
      divBasis = dBasis(:,:,:,1,1); % nExnBxnP
      for d = 2:nC
        divBasis = divBasis + dBasis(:,:,:,d,d); % nExnBxnP
      end
      R = obj.integrate(divBasis, divBasis, k);
    end
    function R = getScaling(obj, nRef)
      assert(strcmp(obj.fesTrial.element.conformity, 'HDiv'));
      R = 2^((nRef*obj.fesTrial.element.dimension));
    end
  end
end
