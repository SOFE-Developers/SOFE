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
      basisI = obj.fesTest.evalGlobalBasis([], obj.codim, 0, {k});
      R = obj.integrate(true, basisI, basisJ, k);
    end
    function R = getScaling(obj, nRef)
      switch obj.fesTrial.element.conformity
        case {'H1', 'L2'}
          R = 2^((nRef*(0-obj.fesTrial.element.dimension)));
        case 'HDiv'
          R = 2^((nRef*(obj.fesTrial.element.dimension-2)));
        case 'HRot'
          R = 2^((nRef*(2-obj.fesTrial.element.dimension)));
      end
    end
  end
end