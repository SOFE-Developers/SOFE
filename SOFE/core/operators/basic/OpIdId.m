classdef OpIdId < Operator % ( c*U, V )
  methods % constructor
    function obj = OpIdId(coeff, codim, fesTrial, varargin) % [fesTest loc]
      obj = obj@Operator(coeff, fesTrial, varargin{:});
      obj.codim = codim;
      if codim == 1 && nargin < 6
        obj.loc = @(x)~obj.fesTest.fixB(x);
      end
    end
  end
  methods
    function R = assembleOp(obj, k)
      points = obj.fesTrial.getQuadData(obj.codim);
      if isempty(points)
        R = obj.fesTrial.mesh.evalFunction(obj.data, points, obj.state, {k}); % nExnP
      else
        basisJ = obj.fesTrial.evalGlobalBasis([], obj.codim, 0, {k}); % nExnBxnPxnC
        basisI = obj.fesTest.evalGlobalBasis([], obj.codim, 0, {k});
        R = obj.integrate(true, basisI, basisJ, k);
      end
    end
  end
end
