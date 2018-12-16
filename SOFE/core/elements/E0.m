classdef E0 < Element
  methods % constructor
    function obj = E0()
      obj = obj@Element(0, 1, 1, 0);
      obj.doFTuple = 1;
      obj.conformity = 'L2';
      obj.isLagrange = true;
    end
  end
  methods % local evaluation
    function B = evalD0Basis(obj, points)
      B = 1;
    end
    function B = evalD1Basis(obj, points)
      B = 0;
    end
  end
end