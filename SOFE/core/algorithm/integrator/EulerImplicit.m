classdef EulerImplicit < Integrator
  methods % constructor
    function obj = EulerImplicit(M0, A, timeline, varargin) % [initCond]
      obj = obj@Integrator(M0, A, timeline, varargin{:});
    end
  end
  methods % integrate
    function compute(obj)
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
        S = S(obj.A.fTest, obj.A.fTrial);
      else
        S = @(x) obj.M0.applySystem(x, obj.A.fTest, obj.A.fTrial) + ...
                 obj.dt(k)*obj.A.applySystem(x, 1);
        b = obj.dt(k)*obj.A.loadVec + obj.M0.applySystem(obj.history{k});
        b = b - obj.M0.applySystem(obj.A.shift) - obj.dt(k)*obj.A.applySystem(obj.A.shift); 
      end
      b = b(obj.A.fTest);      
      obj.history{k+1} = obj.A.shift;
      if obj.directSolve
        if obj.directSolve==2
          if isempty(obj.cache)
            [obj.cache.L,obj.cache.U,obj.cache.P,obj.cache.Q,obj.cache.R] = lu(S);
          end
          obj.history{k+1}(obj.A.fTrial) = obj.cache.Q*(obj.cache.U\(obj.cache.L\(obj.cache.P*(obj.cache.R\b))));
        else
          obj.history{k+1}(obj.A.fTrial) = S\b;
        end
      else
        obj.history{k+1}(obj.A.fTrial) = bicgstab(S, b, 1e-10, 1000, [], [], obj.history{k}(obj.A.fTrial));
      end
    end
  end
end