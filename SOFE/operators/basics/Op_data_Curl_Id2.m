classdef Op_data_Curl_Id2 < Operator % ( c*curl(U), V )
  methods % constructor
    function obj = Op_data_Curl_Id2(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)      
      points = obj.feSpaceTrial.getQuadData(0);
      I = obj.feSpaceTrial.mesh.getBlock(0, k);
      curlBasisJ = obj.feSpaceTrial.evalGlobalBasis(points, [], 'curl', I(1):I(2)); % nExnBxnPxnCxnD
      basisI = obj.feSpaceTest.evalGlobalBasis([], 0, 0, k); % nExnBxnPxnC
      R = obj.integrate(true, basisI, curlBasisJ, k);
    end
  end
end
