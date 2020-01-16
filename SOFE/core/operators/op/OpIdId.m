classdef OpIdId < Operator % ( c*U, V )
  methods % constructor
    function obj = OpIdId(coeff, codim, fesTrial, varargin) % [fesTest loc]
      obj = obj@Operator(coeff, fesTrial, varargin{:});
      obj.codim = codim;
      if codim == 1 && isempty(obj.loc)
        obj.loc = @(x)~obj.fesTest.fixB(x);
      end
    end
  end
  methods
    function R = assembleOp(obj, k)
      basisJ = obj.fesTrial.evalGlobalBasis([], obj.codim, 0, {k}); % nExnBxnPxnC
      if obj.fesTrial ~= obj.fesTest
        basisI = obj.fesTest.evalGlobalBasis([], obj.codim, 0, {k});
        R = obj.integrate(basisI, basisJ, k);
      else
        R = obj.integrate(basisJ, basisJ, k);
      end
    end
    function R = getScaling(obj, nRef)
      dim = obj.fesTrial.element.dimension;
      switch obj.fesTrial.element.conformity
        case {'H1', 'L2'}
          I = 0;
        case 'HDiv'
          I = dim - 1;
        case 'HRot'
          I = 1;
      end
      switch obj.fesTest.element.conformity
        case {'H1', 'L2'}
          J = 0;
        case 'HDiv'
          J = dim - 1;
        case 'HRot'
          J = 1;
      end
      R = 2^(nRef*(I+J-dim));
    end
  end
end