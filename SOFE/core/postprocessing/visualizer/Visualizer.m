classdef Visualizer < SOFE
  properties
    feSpace
  end
  methods % constructor & more
    function obj = Visualizer(feSpace)
      obj.feSpace = feSpace;
    end
    function test(obj, U)
      if numel(U)~=obj.feSpace.getNDoF
        error('w is no DoFVector!');
      end
    end
  end
  methods % visualisation
    function h = show(obj, U, varargin) % [type]
      if nargin > 2
        type = varargin{1};
        if nargin > 3, opt = varargin(2); else, opt = {}; end
        switch type
          case {'g', 'G', 'global', 'surf', 1}
            h = obj.surf(U, opt{:});
          case {'l', 'L', 'local', 'scatter', 0}
            h = obj.scatter(U, opt{:});
          case {'fh', 'FH'}
            h = obj.surfFH(U, opt{:});
          otherwise
            h = obj.patch(U, opt{:});
        end
      else
        h = obj.patch(U);
      end
    end
  end
  methods(Static = true)
    function R = create(fes)
      switch fes.element.dimension
        case 1
          R = Visualizer1D(fes);
        case 2
          R = Visualizer2D(fes);
        case 3
          R = Visualizer3D(fes);
      end
    end
  end
end