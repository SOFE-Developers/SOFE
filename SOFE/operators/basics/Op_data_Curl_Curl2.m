classdef Op_data_Curl_Curl2 < Operator % ( c*curl(U), curl(V) )
  methods % constructor
    function obj = Op_data_Curl_Curl2(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      points = obj.feSpaceTrial.getQuadData(0);
      I = obj.feSpaceTrial.mesh.getBlock(0, k);
      keyboard
      curlBasisJ = obj.feSpaceTrial.evalGlobalBasis(points, [], 'curl', I(1):I(2)); % nExnBxnPxnCxnD
      if obj.isGalerkin
        curlBasisI = curlBasisJ;
      else
        curlBasisI = obj.feSpaceTest.evalGlobalBasis(points, [], 'curl', I(1):I(2)); % nExnBxnPxnCxnD
      end
      R = obj.integrate(true, curlBasisI, curlBasisJ, k);
    end
  end
end
