classdef Static < StaticAlg
  methods % constructor
    function obj = Static(pde, varargin) % [solver]
      obj = obj@StaticAlg(pde, varargin{:});
    end
  end
end
