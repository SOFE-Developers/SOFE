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
        obj.A.setState(obj.timeline.nodes(k+1), obj.history{k+1}); 
        obj.A.assemble();
        [freeI, freeJ] = obj.A.getFreeDoFs;
        b = obj.dt(k)^2*(obj.A.loadVec - obj.A.stiffMat*obj.history{k}) + ...
            M*(2*obj.history{k} - obj.history{k-1}) - ...
            M*obj.A.getShift();
        % solve
        obj.history{k+1} = zeros(size(obj.history{k}));
        obj.history{k+1}(~freeJ) = obj.A.getShift()(~freeJ);
        obj.history{k+1}(freeJ) = obj.A.getShift()(freeJ) + M(freeI, freeJ)\b(freeI);
        if obj.directSolve
          obj.history{k+1}(freeJ) = obj.A.getShift()(freeJ) + M(freeI, freeJ)\b(freeI);
        else
          obj.history{k+1}(freeJ) = obj.A.getShift()(freeJ) + bicgstab(M(freeI, freeJ), b(freeI), 1e-10, 1000, [], [], obj.history{k}(freeJ));
        end
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(tt));
      end
    end
  end
end