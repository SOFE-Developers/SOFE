classdef EulerImplicit2 < Integrator2
  methods % constructor
    function obj = EulerImplicit2(M0, A, timeline, varargin) % [initCond]
      obj = obj@Integrator2(M0, A, timeline, varargin{:});
    end
  end
  methods % integrate
    function integrate(obj, varargin)
      % initial condition
      obj.history{1} = cell(obj.A.nEq,1);
      for k = 1:obj.A.nEq
        try
          obj.history{1}{k} = obj.A.fesTrial{k}.getInterpolation(obj.initCond{k}, 0);
        catch
          obj.history{1}{k} = zeros(obj.A.J{k}(2)-obj.A.J{k}(1)+1, 1);
        end
      end
      obj.history{1} = cell2mat(obj.history{1});
      %
      obj.M0.assemble();
      for k = 1:obj.nT-1
        tt = tic;
        obj.A.setTime(obj.timeline.nodes(k+1));
        obj.A.setState(obj.history{k});
        obj.A.assemble();
        obj.solve(k);
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(tt));
      end
    end
    function solve(obj, k)
      % solve timestep
      if nnz(obj.A.stiffMat)>0 
        S = obj.M0.stiffMat + obj.dt(k)*obj.A.stiffMat;
        b = obj.dt(k)*obj.A.loadVec + obj.M0.stiffMat*obj.history{k};
        b = b - S*obj.A.shift;
        S = S(obj.A.fDoFsTest, obj.A.fDoFsTrial);
      else
        S = @(x) obj.M0.applySystem(x, obj.A.fDoFsTest, obj.A.fDoFsTrial) + ...
                 obj.dt(k)*obj.A.applySystem(x, 1);
        b = obj.dt(k)*obj.A.loadVec + obj.M0.applySystem(obj.history{k});
        b = b - obj.M0.applySystem(obj.A.shift) - obj.dt(k)*obj.A.applySystem(obj.A.shift); 
      end
      b = b(obj.A.fDoFsTest);      
      obj.history{k+1} = obj.A.shift;
      obj.history{k+1}(obj.A.fDoFsTrial) = bicgstab(S, b, 1e-10, 1000, [], [], obj.history{k}(obj.A.fDoFsTrial));
    end
  end
end