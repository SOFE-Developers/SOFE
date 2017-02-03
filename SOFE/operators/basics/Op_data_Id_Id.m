classdef Op_data_Id_Id < Operator % ( c*U, V )
  methods % constructor
    function obj = Op_data_Id_Id(coeff, codim, feSpaceTrial, varargin) % [feSpaceTest loc]
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
      obj.codim = codim;
      if codim == 1 && nargin < 6
        obj.loc = @(x)~obj.feSpaceTest.fixB(x);
      end
    end
  end
  methods
    function R = assembleOp(obj, k)
      points = obj.feSpaceTrial.getQuadData(obj.codim);
      if isempty(points)
        R = obj.feSpaceTrial.mesh.evalFunction(obj.data, points, obj.state, k); % nExnP
      else
        basisJ = obj.feSpaceTrial.evalGlobalBasis([], obj.codim, 0, k); % nExnBxnPxnC
        if obj.isGalerkin
          basisI = basisJ;
        else
          basisI = obj.feSpaceTest.evalGlobalBasis([], obj.codim, 0, k);
        end
        R = obj.integrate(true, basisI, basisJ, k);
      end
    end
  end
end
