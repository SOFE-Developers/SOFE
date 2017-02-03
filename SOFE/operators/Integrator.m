classdef Integrator < SOFEClass
  properties
    massOp, statOp
    mesh
    initCond, nT, dt
    solution
    solver
  end
  methods % constructor
    function obj = Integrator(massOp, statOp, mesh, varargin)
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
end