classdef OpNull < Operator % ( 0*U, V )
  methods 
    function obj = OpNull(fesTrial, fesTest)
      obj = obj@Operator([], fesTrial, fesTest);
    end
  end
  methods
    function R = assembleOp(obj, k)
      nE = obj.fesTest.mesh.getBlock(0, k);
      nD = obj.fesTest.element.dimension;
      nBI = obj.fesTest.element.nB(nD);
      nBJ = obj.fesTrial.element.nB(nD);
      R = zeros(nE(2)-nE(1)+1, nBI, nBJ);
    end
  end
end
