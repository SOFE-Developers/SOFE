classdef OpIdDiv < Operator % (c*u, div(V) )
  methods % constructor
    function obj = OpIdDiv(coeff, fesTrial, fesTest)
      obj = obj@Operator(coeff, fesTrial, fesTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      basisJ = obj.fesTrial.evalGlobalBasis([], 0, 0, {k}); % nExnBxnP
      dBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      nD = size(dBasisI,5);
      divBasisI = dBasisI(:,:,:,1,1); % nExnBxnP
      for d = 2:nD
        divBasisI = divBasisI + dBasisI(:,:,:,d,d); % nExnBxnP
      end
      R = obj.integrate(divBasisI, basisJ, k);
    end
    function R = getScaling(obj, nRef)
      dim = obj.fesTrial.element.dimension;
      switch obj.fesTrial.element.conformity
        case {'H1', 'L2'}
          I = 0;
        otherwise
          assert(0, 'not allowed');
      end
      switch obj.fesTest.element.conformity
        case {'H1', 'L2'}
          J = 1;
        case 'HDiv'
          J = dim;
        otherwise
          assert(0, 'not allowed');
      end
      R = 2^(nRef*(I+J-dim));
    end
  end
end
