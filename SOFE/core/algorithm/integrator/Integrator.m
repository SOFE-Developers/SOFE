classdef Integrator < Algorithm
  properties
    tStep
    timeline
    initCond, nT, dt
    history
  end
  methods % constructor
    function obj = Integrator(timeline, tStep, varargin) % [initCond]
      obj = obj@Algorithm(tStep.pde);
      obj.timeline = timeline;
      obj.dt = diff(timeline.nodes);
      obj.nT = numel(timeline.nodes);
      obj.tStep = tStep;
      if ~isempty(varargin)
        obj.initCond = varargin{1};
        if ~isempty(obj.initCond) && ~iscell(obj.initCond)
          obj.initCond = {obj.initCond};
        end
      end
      obj.history = cell(1,obj.nT);
    end
  end
  methods % integration
    function compute(obj)
      obj.tStep.setFreeDoFs();
      A = obj.tStep.pde;
      % initial condition
      obj.history{1} = cell(A.nEq,1);
      for k = 1:A.nEq
        if ~isempty(obj.initCond)
          obj.history{1}{k} = A.fesTrial{k}.getInterpolation(obj.initCond{k}, 0);
        else
          obj.history{1}{k} = zeros(A.J(k,2)-A.J(k,1)+1, 1);
        end
      end
      obj.history{1} = cell2mat(obj.history{1});
      % starting values
      kS = 1;
      if obj.tStep.nS>1
        for kS = 2:obj.tStep.nS
          obj.history{kS} = obj.history{1};
        end
      end
      % time loop
      for k = kS:obj.nT-1
        cTime = tic;
        obj.history{k+1} = obj.tStep.compute(obj.timeline.nodes([k, k+1]), cell2mat(obj.history(k-obj.tStep.nS+1:k)));
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(cTime));
      end
    end
    function compute2(obj)
      % routine for space adaptivity 
      A = obj.tStep.pde;
      % initial condition
      obj.history{1} = cell(A.nEq,1);
      for k = 1:A.nEq
        if ~isempty(obj.initCond)
          obj.history{1}{k} = A.fesTrial{k}.getInterpolation(obj.initCond{k}, 0);
        else
          obj.history{1}{k} = zeros(A.J(k,2)-A.J(k,1)+1, 1);
        end
      end
      obj.history{1} = cell2mat(obj.history{1});
      assert(obj.tStep.nS==1);
      v = Visualizer.create(obj.pde.fesTrial{1});
      % time loop
      for k = 1:obj.nT-1
        cTime = tic;
        %
        obj.tStep.M0.assemble();
        obj.tStep.setFreeDoFs();
        %
        obj.history{k+1} = obj.tStep.compute(obj.timeline.nodes([k, k+1]), obj.history{k});
        %
%         figure(1),clf, v.show(obj.history{k+1},'p'); drawnow
%         figure(2),obj.pde.mesh.show(); drawnow
        %
        for cc = 1:3
%           I = obj.pde.fesTrial{1}.evalDoFVector(obj.history{k+1},[1 1]/3,[],0)>0.02;
          I =  max(obj.history{k+1}(obj.pde.mesh.topology.getEntity('0')),[],2)>0.02;
          I = I & obj.pde.mesh.getMeasure('0')>0.005^2+1e-12;
          if isempty(I), break; end
          obj.history{k+1} = obj.pde.mesh.adaptiveRefine(I)*obj.history{k+1};
        end
        for cc = 1:3
%           I = obj.pde.fesTrial{1}.evalDoFVector(obj.history{k+1},[1 1]/3,[],0)<0.02;
          I =  max(obj.history{k+1}(obj.pde.mesh.topology.getEntity('0')),[],2)<0.02;
          I = I & obj.pde.mesh.getMeasure('0')<0.0125^2+1e-12;
          if isempty(I), break; end
          obj.history{k+1} = obj.pde.mesh.coarsen(I)*obj.history{k+1};  
        end
        %
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(cTime));
      end
    end
  end
  methods(Static=true) % Test
    function R = testConvergence(intType, dim, pSpace)
      doVis = 0;
      T = 1; epsilon = 1e-0; omega = 1; isHom = 0;
      switch dim
        case 1
          if isHom
            uEx = @(x,t)sin(pi*x(:,1)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*sin(pi*x(:,1)).*cos(2*pi*omega*t) + (pi)^2*epsilon*uEx(x,t);
          else
            uEx = @(x,t)cos(pi/2*x(:,1)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*cos(pi/2*x(:,1)).*cos(2*pi*omega*t) + pi^2/4*epsilon*uEx(x,t);
          end
          qr = GaussInt(5);
        case 2
          if isHom
            uEx = @(x,t)sin(pi*x(:,1)).*sin(pi*x(:,2)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*sin(pi*x(:,1)).*sin(pi*x(:,2)).*cos(2*pi*omega*t) + 2*pi^2*epsilon*uEx(x,t);
          else
            uEx = @(x,t)cos(pi/2*x(:,1)).*sin(pi*x(:,2)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*cos(pi/2*x(:,1)).*sin(pi*x(:,2)).*cos(2*pi*omega*t) + (pi^2/4 + pi^2)*epsilon*uEx(x,t);
          end
         %
         qr = GaussQuad(5);
      end
%       NN = 50;
      NN = (10:20:70);
%       NN = [5 10 15 20];
      R = zeros(numel(NN),1); cnt = 1;
      for N = NN
        M = N;
        % SOLVE
        m = RegularMesh(repmat(N,dim,1), repmat([0 1],dim,1), 0);
        fes = FESpace(m, QpL(dim, pSpace), @(x)x(:,1)<Inf, uEx);
        p = Poisson(struct('a',epsilon, 'f',f), fes);
        q = Integrator(RegularMesh(M, [0 T]), TimeStep.create(intType, Mass(fes), p, DirectSol()));
        q.compute();
        % COMPUTE ERROR
        err = zeros(M+1,1);
        for k = 1:M+1
          UEX = m.evalFunction(@(x)uEx(x,T*(k-1)/M),qr.points,[]);
          UU = fes.evalDoFVector(q.history{k}, qr.points,[],0);
          err(k) = abs(m.integrate((UEX - UU).^2, qr)).^0.5; % L2-norm
%           err(k) = max(max(abs(UEX - UU))); % max-norm
        end
        R(cnt) = max(err);
        cnt = cnt + 1;
      end
      clf
      loglog(NN,R); hold on
      loglog(NN,1./NN.^1);
      loglog(NN,1./NN.^2);
      loglog(NN,1./NN.^3);
      loglog(NN,1./NN.^4);
      loglog(NN,1./NN.^5);
      loglog(NN,1./NN.^6); hold off
      %
      if ~doVis, return; end
      v = Visualizer.create(fes);
      for k = 1:1:q.nT
        clf
        switch dim
          case 1
            v.show(q.history{k}, 'p');
            axis([0 1 -1 1])
          case 2
            v.show(q.history{k}, 'p');
            axis([0 1 0 1 -1 1])
            caxis([-1 1]); 
            view(3);
        end
        fprintf('timestep: %d / %d\n', k, q.nT);
        title(num2str((k-1)*T/M))
        pause(0.01)
      end
    end
  end
end