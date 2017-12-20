classdef Integrator2 < SOFE
  properties
    M0, A
    timeline
    initCond, nT, dt
    history
  end
  methods % constructor
    function obj = Integrator2(M0, A, timeline, varargin) % [initCond]
      obj.M0 = M0; obj.A = A; obj.timeline = timeline;
      obj.dt = diff(timeline.nodes); obj.nT = numel(timeline.nodes);
      if nargin > 3
        obj.initCond = varargin{1};
        if ~iscell(obj.initCond), obj.initCond = {obj.initCond}; end
      end
      obj.history = cell(obj.nT,1);
    end
  end
  methods % integrate
    function integrate(obj) %#ok<MANU>
    end
  end
end