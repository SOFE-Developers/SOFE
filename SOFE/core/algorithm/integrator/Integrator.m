classdef Integrator < Algorithm
  properties
    tStep
    timeline
    initCond, nT, dt
    history
    fullHistory
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
      obj.fullHistory = cell(1,obj.nT);
    end
  end
  methods % integration
    function compute(obj)
      obj.tStep.setFreeDoFs();
      A = obj.tStep.pde;
      % initial condition
      obj.history{1} = cell(A.nEq,1);
      if ~iscell(obj.initCond), obj.initCond = {obj.initCond}; end
      if numel(obj.initCond) < A.nEq, obj.initCond = repmat(obj.initCond, 1,A.nEq); end
      for k = 1:A.nEq
        if ~isempty(obj.initCond{k})
          obj.history{1}{k} = A.fesTrial{k}.getL2Interpolant(obj.initCond{k}, A.fesTrial{k}.mesh.element.dimension);
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
        obj.fullHistory{k+1} = obj.tStep.compute(obj.timeline.nodes([k, k+1]), cell2mat(obj.history(k-obj.tStep.nS+1:k))); % check cGP!!!
        try
          obj.history{k+1} = obj.tStep.eval(obj.fullHistory{k+1}, 1.0);
        catch
          obj.history{k+1} = obj.fullHistory{k+1}(:,1);
        end
        fprintf('timestep: %d / %d: %f sec\n', k, obj.nT-1, toc(cTime));
      end
      % constant "approximation" of initial value
      obj.fullHistory{1} = zeros(size(obj.fullHistory{2})); 
      obj.fullHistory{1}(:,1) = obj.history{1}; % left linear part
      obj.fullHistory{1}(:,2) = obj.history{1}; % right linear part
    end
  end
  methods % access
    function sol = getSolution(obj,varargin)
      if nargin>1
        nr = varargin{1};
      else
        nr = 1;
      end
      sol = cell(1,numel(obj.fullHistory));
      if isnumeric(nr)
        id = 0;
        for k=1:nr-1
          id = id+obj.pde.fesTrial{k}.getNDoF;
        end
        id = id+(1:obj.pde.fesTrial{nr}.getNDoF);
      else
        id = ':';
      end
      for i=1:numel(obj.fullHistory)
        sol{i} = obj.fullHistory{i}(id,:);
      end
    end
    function sol = setSolution(obj,sol,nr)
      try
        if numel(obj.fullHistory)~=numel(sol)
          obj.fullHistory = cell(1,numel(sol));
          obj.history     = cell(1,numel(sol));
        end
      catch
        obj.fullHistory = cell(1,numel(sol));
        obj.history     = cell(1,numel(sol));
      end
      if isnumeric(nr)
        id = 0;
        for k=1:nr-1
          id = id+obj.pde.fesTrial{k}.getNDoF;
        end
        id = id+(1:obj.pde.fesTrial{nr}.getNDoF);
      else
        id = ':';
      end
      for i=1:numel(sol)
        obj.fullHistory{i}(id,:)=sol{i};
        obj.history{i}(id,:)=sol{i}(:,end);
      end
    end
  end
  methods(Static=true) % Test
    function R = testConvergence(intType, dim, pSpace, varargin)
      try doVis = varargin{1}.doVis; catch, doVis = false; end
      try T = varargin{1}.T; catch, T = 1; end
      try NVec = varargin{1}.NVec(:); catch, NVec = (10:10:50)'; end
      try omega = varargin{1}.omega; catch, omega = 1; end
      try isHom = varargin{1}.isHom; catch, isHom = true; end
      try aTime = varargin{1}.aTime; catch, aTime = false; end
      try mTime = varargin{1}.mTime; catch, mTime = false; end
      if aTime
        aCoeff = @(x,t)(1+0.5*sin(2*pi/T*t)) + 0*x(:,1);
      else
        aCoeff = @(x,t)1+0*x(:,1);
      end
      if mTime
        mCoeff = @(x,t)(1+0.5*cos(2*pi/T*t)) + 0*x(:,1);
      else
        mCoeff = @(x,t)1+0*x(:,1);
      end
      switch dim
        case 1
          if isHom
            uEx = @(x,t)sin(pi*x(:,1)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*sin(pi*x(:,1)).*cos(2*pi*omega*t).*mCoeff(x,t) + ...
                      (pi)^2*aCoeff(x,t).*uEx(x,t);
          else
            uEx = @(x,t)cos(pi/2*x(:,1)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*cos(pi/2*x(:,1)).*cos(2*pi*omega*t).*mCoeff(x,t) + ...
                      pi^2/4*aCoeff(x,t).*uEx(x,t);
          end
          qr = GaussInt(2*pSpace+1);
        case 2
          if isHom
            uEx = @(x,t)sin(pi*x(:,1)).*sin(pi*x(:,2)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*sin(pi*x(:,1)).*sin(pi*x(:,2)).*cos(2*pi*omega*t).*mCoeff(x,t) + ...
                      2*pi^2*aCoeff(x,t).*uEx(x,t);
          else
            uEx = @(x,t)cos(pi/2*x(:,1)).*sin(pi*x(:,2)).*sin(2*pi*omega*t);
            f = @(x,t)2*pi*omega*cos(pi/2*x(:,1)).*sin(pi*x(:,2)).*cos(2*pi*omega*t).*mCoeff(x,t) + ...
                      (pi^2/4 + pi^2)*aCoeff(x,t).*uEx(x,t);
          end
          qr = GaussQuad(2*pSpace+1);
      end
      R = zeros(numel(NVec),1); cnt = 1;
      for N = NVec'
        M = N;
        % SOLVE
        m = RegularMesh(repmat(N,dim,1), repmat([0 1],dim,1), 0);
        fes = FESpace(m, QpL(dim, pSpace), @(x)x(:,1)<Inf, uEx);
        p = Poisson(struct('a',aCoeff, 'f',f), fes);
        q = Integrator(RegularMesh(M, [0 T]), TimeStep.create(intType, Mass2({mCoeff}, {fes}), p, DirectSol()));
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
      loglog(NVec,R); 
      hold on
      for k=1:7
        loglog(NVec,1./NVec.^k,'k');
      end
      hold off
      R = [R [0; -log2(R(2:end)./R(1:end-1))./log2(NVec(2:end)./NVec(1:end-1))]];
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