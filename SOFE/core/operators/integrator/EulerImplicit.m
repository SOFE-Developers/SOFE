classdef EulerImplicit < Integrator
  methods % constructor
    function obj = EulerImplicit(massOp, statOp, mesh, varargin) % [initCond]
      obj = obj@Integrator(massOp, statOp, mesh, varargin{:});
    end
  end
  methods % integrate
    function integrate(obj, varargin)
      % initial condition
      obj.solution{1} = cell(obj.statOp.nEq,1);
      for k = 1:obj.statOp.nEq
        try
          obj.solution{1}{k} = obj.statOp.fesTrial{k}.getInterpolation(obj.initCond{k}, 0);
        catch
          obj.solution{1}{k} = zeros(obj.statOp.J{k}(2)-obj.statOp.J{k}(1)+1, 1);
        end
      end
      obj.solution{1} = cell2mat(obj.solution{1});
      %
      obj.massOp.assemble();
      for k = 1:obj.nT-1
        tt = tic;
        obj.statOp.setTime(obj.mesh.nodes(k+1));
        obj.statOp.setState(obj.solution{k});
        % assemble
        obj.statOp.assemble();
        S = obj.massOp.stiffMat + obj.dt(k)*obj.statOp.stiffMat;
        b = obj.dt(k)*obj.statOp.loadVec + obj.massOp.stiffMat*obj.solution{k} - S*obj.statOp.shift;
        iTest = obj.statOp.fDoFsTest; iTrial = obj.statOp.fDoFsTrial;
        % solve
        obj.solution{k+1} = zeros(size(obj.solution{k}));
        obj.solution{k+1}(~iTrial) = obj.statOp.shift(~iTrial);
        obj.solution{k+1}(iTrial) = obj.statOp.shift(iTrial) + ...
          obj.solver.solve(S(iTest, iTrial), b(iTest), obj.solution{k}(iTrial));
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(tt));
      end
    end
  end
end