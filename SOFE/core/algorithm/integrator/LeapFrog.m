classdef LeapFrog < Integrator
  methods % constructor
    function obj = LeapFrog(M0, A, timeline, varargin) % [initCond]
      obj = obj@Integrator(M0, A, timeline, varargin{:});
    end
  end
  methods % integrate
    function compute(obj)
      % initial condition
      obj.history{1} = cell(obj.A.nEq,1);
      for k = 1:obj.A.nEq
        if ~isempty(obj.initCond)
          obj.history{1}{k} = obj.A.fesTrial{k}.getInterpolation(obj.initCond{k}, 0);
        else
          obj.history{1}{k} = zeros(obj.A.dimTrial(k), 1);
        end
      end
      obj.history{1} = cell2mat(obj.history{1});
      obj.history{2} = obj.history{1};
      %
      obj.M0.assemble();
      M = obj.M0.stiffMat;
      for k = 2:obj.nT-1
        tt = tic;
        obj.A.setTime(obj.timeline.nodes(k+1));
        obj.A.setState(obj.history{k+1}); 
        obj.A.assemble();
        iTest = obj.A.fTest; iTrial = obj.A.fTrial;
        b = obj.dt(k)^2*(obj.A.loadVec - obj.A.stiffMat*obj.history{k}) + ...
            M*(2*obj.history{k} - obj.history{k-1}) - ...
            M*obj.A.shift;
        % solve
        obj.history{k+1} = zeros(size(obj.history{k}));
        obj.history{k+1}(~iTrial) = obj.A.shift(~iTrial);
        obj.history{k+1}(iTrial) = obj.A.shift(iTrial) + M(iTest, iTrial)\b(iTest);
        if obj.directSolve
          obj.history{k+1}(iTrial) = obj.A.shift(iTrial) + M(iTest, iTrial)\b(iTest);
        else
          obj.history{k+1}(iTrial) = obj.A.shift(iTrial) + bicgstab(M(iTest, iTrial), b(iTest), 1e-10, 1000, [], [], obj.history{k}(obj.A.fTrial));
        end
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(tt));
      end
    end
  end
end