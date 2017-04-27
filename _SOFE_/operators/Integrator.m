classdef Integrator < SOFE
  properties
    massOp, statOp
    mesh
    solver
    initCond, nT, dt
    solution
  end
  methods % constructor
    function obj = Integrator(massOp, statOp, mesh, varargin) % [initCond]
      obj.solver = statOp.solver;
      obj.massOp = massOp;
      obj.statOp = statOp;
      obj.mesh = mesh;
      obj.dt = diff(mesh.topology.nodes);
      obj.nT = numel(mesh.topology.nodes);
      if nargin > 3
        obj.initCond = varargin{1};
        if ~iscell(obj.initCond), obj.initCond = {obj.initCond}; end
      end
      obj.solution = cell(obj.nT,1);
    end
  end
  methods % integrate
    function integrate(obj)
    end
  end
end