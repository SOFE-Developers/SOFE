classdef Visualizer1D < Visualizer
  methods % constructor
    function obj = Visualizer1D(feSpace)
      obj = obj@Visualizer(feSpace);
    end
  end
  methods % display
    function h = patch(obj, U, varargin)
      obj.test(U);
      X = sort(obj.feSpace.mesh.topology.getEntity(0));
      Y = obj.feSpace.evalDoFVector(U, {X}, [], 0);
      h = plot(X,Y);
    end
    function h = surf(obj, U, varargin)
      obj.test(U);
      try N = varargin{1}.N; catch, N = 200; end
      try
        box = varargin{1}.box;
      catch
        gs = obj.feSpace.mesh.topology.getGlobalSearcher();
        box = gs.diam';
      end
      X = linspace(box(1), box(2), N)';
      Y = obj.feSpace.evalDoFVector(U, {X}, [], 0);
      h = plot(X,Y);
    end
    function h = scatter(obj, U, varargin) % [resolution]
      if numel(U)~=obj.feSpace.getNDoF
        error('w is no DoFVector!');
      end
      if nargin > 2, res = varargin{:}; else, res = 10; end
      points = linspace(0,1,res)';
      P = obj.feSpace.mesh.evalReferenceMap(points, 0); % nExnPxnD
      Z = obj.feSpace.evalDoFVector(U, points, [], 0);
      if size(P,3) == 1
        h = scatter(P(:), Z(:), 2);
      elseif size(P,3) == 2
        Px = P(:,:,1); Py = P(:,:,2);
        h = plot3k([Px(:),Py(:),Z(:)]); view(0,90);
      else
        Px = P(:,:,1); Py = P(:,:,2); Pz = P(:,:,3);
        h = plot3k([Px(:),Py(:), Pz(:)], 'ColorData', Z(:));
      end
    end
  end
end