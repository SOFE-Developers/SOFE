classdef ThetaMethod < Integrator
  properties
    theta = 0.5;
  end
  methods % constructor
    function obj = ThetaMethod(M0, A, timeline, varargin) % [initCond]
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
          obj.history{1}{k} = zeros(obj.A.J(k,2)-obj.A.J(k,1)+1, 1);
        end
      end
      obj.history{1} = cell2mat(obj.history{1});
      %
      obj.M0.assemble();
      obj.A.assemble();
      for k = 1:obj.nT-1
        tt = tic;
        obj.A.setState(obj.timeline.nodes(k+1), obj.history{k});
        loadPrev = obj.A.loadVec;
        obj.A.assemble();
        obj.solve(k, loadPrev);
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(tt));
      end
    end
    function solve(obj, k, loadPrev)
      % solve timestep
      [freeI, freeJ] = obj.A.getFreeDoFs();
      if nnz(obj.A.stiffMat)>0
        if obj.A.narginData>1, warning('F must be constant'); end
        S = (obj.M0.stiffMat+(obj.theta*obj.dt(k))*obj.A.stiffMat);
        b = obj.M0.stiffMat*(obj.history{k}-obj.A.getShift()) + ...
              obj.dt(k)*(obj.theta*(obj.A.loadVec-loadPrev + ...
                         obj.A.stiffMat*(obj.history{k}-obj.A.getShift())) + ...
                         loadPrev - obj.A.stiffMat*obj.history{k});
      else
        error('TODO');
      end
      obj.history{k+1} = obj.A.getShift();
      S = S(freeI, freeJ); b = b(freeI);
      if obj.directSolve
        if obj.directSolve==2
          if isempty(obj.cache)
            [obj.cache.L,obj.cache.U,obj.cache.P,obj.cache.Q,obj.cache.R] = lu(S);
          end
          obj.history{k+1}(freeJ) = obj.cache.Q*(obj.cache.U\(obj.cache.L\(obj.cache.P*(obj.cache.R\b))));
        else
          obj.history{k+1}(freeJ) = S\b;
        end
      else
        obj.history{k+1}(freeJ) = bicgstab(S, b, 1e-10, 1000, [], [], obj.history{k}(freeJ));
      end
    end
  end
end