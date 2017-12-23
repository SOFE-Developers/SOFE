classdef OpIntId < Operator % ( \int(u), lambda )_R
  properties
    subOp
  end
  methods % constructor
    function obj = OpIntId(fes)
      obj = obj@Operator(@(x)1+0*x(:,1), fes);
      obj.subOp = FcId(1, fes, 0);
      obj.fesTest = ScalarFESpace();
    end
  end
  methods
    function assemble(obj)
      obj.subOp.assemble();
      obj.matrix = obj.subOp.matrix';
    end
  end
end
