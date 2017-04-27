classdef Visualizer1D < Visualizer
  methods % constructor
    function obj = Visualizer1D(feSpace)
      obj = obj@Visualizer(feSpace);
    end
  end
  methods % display
    function patch(obj, U, varargin)
      obj.surf(U, varargin{:});
    end
    function surf(obj, U, varargin)
      f = @(x)obj.feSpace.evalDoFVector(U, {x}, [], 0);
      fplot(f, obj.feSpace.mesh.topology.globalSearcher.diam, varargin{:});
    end
    function scatter(obj, U, varargin) % [resolution]
      if numel(U)~=obj.feSpace.getNDoF
        error('w is no DoFVector!');
      end
      if nargin > 3, res = varargin{:}; else res = 10; end
      points = linspace(0,1,res)';
      P = obj.feSpace.mesh.evalReferenceMap(points, 0); % nExnPxnD
      Z = obj.feSpace.evalDoFVector(U, points, [], 0);
      if size(P,3) == 1
        scatter(P(:), Z(:), 2)
      elseif size(P,3) == 2
        Px = P(:,:,1); Py = P(:,:,2);
        plot3k([Px(:),Py(:),Z(:)]); view(0,90);
      else
        Px = P(:,:,1); Py = P(:,:,2); Pz = P(:,:,3);
        plot3k([Px(:),Py(:), Pz(:)], 'ColorData', Z(:));
      end
    end
  end
end