classdef LeapFrog < Integrator
  methods % constructor
    function obj = LeapFrog(massOp, statOp, mesh, varargin) % [initCond]
      obj = obj@Integrator(massOp, statOp, mesh, varargin{:});
    end
  end
  methods % integrate
    function integrate(obj, varargin)
      % initial condition
      obj.solution{1} = cell(obj.statOp.nEq,1);
      for k = 1:obj.statOp.nEq
        if ~isempty(obj.initCond)
          obj.solution{1}{k} = obj.statOp.fesTrial{k}.getInterpolation(obj.initCond{k}, 0);
        else
          obj.solution{1}{k} = zeros(obj.statOp.dimTrial(k), 1);
        end
      end
      obj.solution{1} = cell2mat(obj.solution{1});
      obj.solution{2} = obj.solution{1};
      %
      obj.massOp.assemble();
      M = obj.massOp.stiffMat;
      for k = 2:obj.nT-1
        tt = tic;
        obj.statOp.setTime(obj.mesh.nodes(k+1));
        obj.statOp.assemble();
        iTest = obj.statOp.fDoFsTest; iTrial = obj.statOp.fDoFsTrial;
        b = obj.dt(k)^2*(obj.statOp.loadVec - obj.statOp.stiffMat*obj.solution{k}) + ...
            M*(2*obj.solution{k} - obj.solution{k-1}) - ...
            M*obj.statOp.shift;
        % solve
        obj.solution{k+1} = zeros(size(obj.solution{k}));
        obj.solution{k+1}(~iTrial) = obj.statOp.shift(~iTrial);
        obj.solution{k+1}(iTrial) = obj.statOp.shift(iTrial) + ...
          obj.solver.solve(M(iTest, iTrial), b(iTest));
        obj.statOp.setState(obj.solution{k+1});    
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(tt));
      end
    end
  end
end