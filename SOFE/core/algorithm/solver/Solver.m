classdef Solver < SOFE
  properties
    pde
    solution % deprecated
  end
  methods % constructor
    function obj = Solver(varargin)
      obj = obj@SOFE();
      if ~isempty(varargin), obj.pde = varargin{1}; end
    end
  end
end
