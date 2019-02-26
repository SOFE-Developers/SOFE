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
      colormap hot
    end
  end
  methods(Static = true)
    function R = create(fes)
      R = [];
      dim = fes.mesh.topology.dimP;
      eval(['R = Visualizer' num2str(dim) 'D(fes);']);
    end
    function rotate()
      for k = 1:360
        view(k,ceil(k/4));
        drawnow;
      end
    end
  end
end