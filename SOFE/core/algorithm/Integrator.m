classdef Integrator < Algorithm
  properties
    M0, A
    timeline
    initCond, nT, dt
    history
    %
    directSolve = 0; % 0: off, 1: on, 2: use cache
    cache
  end
  methods % constructor
    function obj = Integrator(M0, A, timeline, varargin) % [initCond]
      obj = obj@Algorithm(A);
      obj.M0 = M0; obj.A = A;
      obj.timeline = timeline;
      obj.dt = diff(timeline.nodes);
      obj.nT = numel(timeline.nodes);
      if nargin > 3
        obj.initCond = varargin{1};
        if ~iscell(obj.initCond), obj.initCond = {obj.initCond}; end
      end
      obj.history = cell(obj.nT,1);
    end
  end
end