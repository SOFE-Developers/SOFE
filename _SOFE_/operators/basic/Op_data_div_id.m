classdef Op_data_div_id < Operator % (c*div(U), v )
  methods % constructor
    function obj = Op_data_div_id(coeff, fesTrial, fesTest)
      obj = obj@Operator(coeff, fesTrial, fesTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      nD = size(dBasisJ,5);
      basisI = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnP
      divBasisJ = dBasisJ(:,:,:,1,1); % nExnBxnP
      for d = 2:nD
        divBasisJ = divBasisJ + dBasisJ(:,:,:,d,d); % nExnBxnP
      end
      R = obj.integrate(true, basisI, divBasisJ, k);
    end
  end
end
