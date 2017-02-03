classdef EulerImplicit < Integrator
  methods % constructor
    function obj = EulerImplicit(massOp, statOp, mesh, varargin)
      obj = obj@Integrator(massOp, statOp, mesh, varargin{:});
    end
  end
  methods
    function integrate(obj, varargin)
      % initial condition
      obj.solution{1} = cell(obj.statOp.nEq,1);
      for k = 1:obj.statOp.nEq
        if ~isempty(obj.initCond)
          obj.solution{1}{k} = obj.statOp.fesTrial{k}.getWeakInterpolation(obj.initCond{k}, 0);
        else
          obj.solution{1}{k} = zeros(obj.statOp.dimTrial(k), 1);
        end
      end
      obj.solution{1} = cell2mat(obj.solution{1});
      %
      obj.massOp.assemble();
      for k = 1:obj.nT-1
        tt = tic;
        obj.statOp.setTime(obj.mesh.topology.nodes(k+1));
        % assemble
        obj.statOp.assemble();
        S = obj.massOp.stiffMat + obj.dt(k)*obj.statOp.stiffMat;
        b = obj.dt(k)*obj.statOp.loadVec + obj.massOp.stiffMat*obj.solution{k} - S*obj.statOp.shift;
        iTest = obj.statOp.fDoFsTest; iTrial = obj.statOp.fDoFsTrial;
        % solve
        obj.solution{k+1} = zeros(size(obj.solution{k}));
        obj.solution{k+1}(~iTrial) = obj.statOp.shift(~iTrial);
        obj.solution{k+1}(iTrial) = obj.statOp.shift(iTrial) + ...
          obj.solver.solve(S(iTest, iTrial), b(iTest));
        obj.statOp.setState(obj.solution{k+1});
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(tt));
      end
    end
  end
end