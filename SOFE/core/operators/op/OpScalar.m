classdef OpScalar < Operator % ( u, v )_R
  properties
    scalar
  end
  methods % constructor
    function obj = OpScalar(scalar)
      obj = obj@Operator(1, ScalarFESpace());
      obj.scalar = scalar;
    end
  end
  methods
    function assemble(obj)
      obj.matrix = obj.scalar;
    end
  end
end
