classdef ScalarFESpace < SOFE
  properties
    narginShift = 0;
  end
  methods
    function obj = ScalarFESpace()
      obj = obj@SOFE();
    end
    function R = getNDoF(obj);
      R = 1;
    end
    function R = getFreeDoFs(obj);
      R = true;
    end
    function R = getShift(obj, varargin)
      R = 0;
    end
    function R = evalDoFVector(varargin)
      R = [];
    end
  end
end