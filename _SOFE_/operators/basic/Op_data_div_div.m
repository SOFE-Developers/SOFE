classdef Op_data_div_div < Operator % ( c*div(U), div(V) )
  methods % constructor
    function obj = Op_data_div_div(coeff, fes)
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
      R = obj.integrate(true, divBasis, divBasis, k);
    end
  end
end
