classdef Solver < SOFE
  properties
    pde
    solution % deprecated
  end
  methods % constructor
    function obj = Solver(varargin)
      obj = obj@SOFE();
      if ~isempty(varargin),
        warning('deprecated');
        obj.pde = varargin{1};
      end
    end
    function setPDE(obj, pde)
      obj.pde = pde;  
    end
  end
end
