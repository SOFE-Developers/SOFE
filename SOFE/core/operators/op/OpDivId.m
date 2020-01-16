classdef OpDivId < Operator % (c*div(U), v )
  methods % constructor
    function obj = OpDivId(coeff, fesTrial, fesTest)
      obj = obj@Operator(coeff, fesTrial, fesTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
%       dBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
%       nD = size(dBasisJ,5);
%       divBasisJ = dBasisJ(:,:,:,1,1); % nExnBxnP
%       for d = 2:nD
%         divBasisJ = divBasisJ + dBasisJ(:,:,:,d,d); % nExnBxnP
%       end
      basisI = obj.fesTest.evalGlobalBasis([], 0, 0, {k}); % nExnBxnP
      divBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 'div', {k}); % nExnBxnP
      R = obj.integrate(basisI, divBasisJ, k);
    end
    function R = getScaling(obj, nRef)
      dim = obj.fesTrial.element.dimension;
      switch obj.fesTrial.element.conformity
        case {'H1', 'L2'}
          I = 1;
        case 'HDiv'
          I = dim;
        otherwise
          assert(0, 'not allowed');
      end
      switch obj.fesTest.element.conformity
        case {'H1', 'L2'}
          J = 0;
        otherwise
          assert(0, 'not allowed');
      end
      R = 2^(nRef*(I+J-dim));
    end
  end
end
