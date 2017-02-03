classdef Op_NULL < Operator % (( 0*U, V ))
  methods % (( 0*U, V ))
    function obj = Op_NULL(feSpaceTrial, feSpaceTest)
      obj = obj@Operator([], feSpaceTrial, feSpaceTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      nE = obj.feSpaceTest.mesh.getBlock(0, k);
      nD = obj.feSpaceTest.element.dimension;
      nBI = obj.feSpaceTest.element.nB(nD);
      nBJ = obj.feSpaceTrial.element.nB(nD);
      R = zeros(nE(2)-nE(1)+1, nBI, nBJ);
    end
  end
end
