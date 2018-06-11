classdef TimeStep < Algorithm
  properties
    data
    M0
    nS
  end
  methods % constructor
    function obj = TimeStep(data, M0, pde, solver)
      obj = obj@Algorithm(pde);
      obj.data = data;
      obj.M0 = M0;
      obj.M0.assemble();
      obj.setSolver(solver);
    end
  end
  methods(Static = true)
    function R = create(intType, M0, pde, solver) % [solver]
      if numel(intType)>4 && strcmp(intType(1:5),'theta')
        th = str2double(intType(6:end));
%        R = RungeKuttaMethod(struct('A',[0 0;1-th th],'b',[1-th th],'c',[0 1]), M0, pde, solver);
        R = Theta(th, M0, pde, solver);
        return
      end
      switch intType
        % explicit Runge Kutta
        case 'EulerEx'
          R = RungeKuttaMethod(struct('A',0,'b',1,'c',0), M0, pde, solver);
        case 'Heun2'
          R = RungeKuttaMethod(struct('A',[0 0;1 0],'b',[1 1]/2,'c',[0 1]), M0, pde, solver);
        case 'Heun3'
          R = RungeKuttaMethod(struct('A',[0 0 0;1/3 0 0; 0 2/3 0],'b',[1 0 3]/4,'c',[0 1 2]/3), M0, pde, solver);
        case 'Ralston'
          R = RungeKuttaMethod(struct('A',[0 0;2/3 0],'b',[1 3]/4,'c',[0 2/3]), M0, pde, solver);
        case 'RK3'
          R = RungeKuttaMethod(struct('A',[0 0 0;0.5 0 0;-1 2 0],'b',[1 4 1]/6,'c',[0 0.5 1]), M0, pde, solver);
        case 'RK4'
          R = RungeKuttaMethod(struct('A',[0 0 0 0;0.5 0 0 0;0 0.5 0 0;0 0 1 0],'b',[1 2 2 1]/6,'c',[0 0.5 0.5 1]), M0, pde, solver);
        case '3/8'
          R = RungeKuttaMethod(struct('A',[0 0 0 0;1 0 0 0;-1 3 0 0;3 -3 3 0]/3,'b',[1 3 3 1]/8,'c',[0 1 2 3]/3), M0, pde, solver);
        case 'RKF'
          R = RungeKuttaMethod(struct('A',[0 0 0 0 0 0;0.25 0 0 0 0 0;3/32 9/32 0 0 0 0;1932/2197 -7200/2197  7296/2197 0 0 0;439/216 -8 3680/513 -845/4104 0 0;-8/27 2 -3544/2565 1859/4104 -11/40 0],'b',[16/135 0 6656/12825 28561/56430 -9/50 2/55],'c',[0 0.25 0.375 12/13 1 0.5]), M0, pde, solver);
        % implicit Runge Kutta
        case 'EulerIm'
%          R = RungeKuttaMethod(struct('A',1,'b',1,'c',1), M0, pde, solver);
          R = EulerImp(M0, pde, solver);
        case 'CN'
          R = RungeKuttaMethod(struct('A',[0 0; 0.5 0.5],'b',[0.5 0.5],'c',[0 1]), M0, pde, solver);
        case 'Gauss2'
          % order 4
          R = RungeKuttaMethod(struct('A',[0.25 0.25-sqrt(3)/6; 0.25+sqrt(3)/6 0.25],'b',[0.5 0.5],'c',[0.5-sqrt(3)/6 0.5+sqrt(3)/6]), M0, pde, solver);
        case 'Gauss3'
          % order 6
          R = RungeKuttaMethod(struct('A',[5/36 2/9-sqrt(15)/15 5/36-sqrt(15)/30; 5/36+sqrt(15)/24 2/9 5/36-sqrt(15)/24; 5/36+sqrt(15)/30 2/9+sqrt(15)/15 5/36],'b',[5/18 4/9 5/18],'c',[0.5-sqrt(15)/10 0.5 0.5+sqrt(15)/10]), M0, pde, solver);
        case 'SDIRK3'
          % order 4
          s0 = 1; % 0/0.5/1
          s = fsolve(@(x)x.^3-1.5*x.^2+0.5*x-1/24,s0);
          b = 1/(6*(1-2*s)^2);
          R = RungeKuttaMethod(struct('A',[s 0 0;0.5-s s 0;2*s 1-4*s s],'b',[b 1-2*b b],'c',[s 0.5 1-s]), M0, pde, solver);
        case 'SDIRK(3,4,5)'
          s = -1.068579021301629;
          b = 0.6696236404609742;
          c = 0.08902038200616;
          d = 0.702767557525405;
          R = RungeKuttaMethod(struct('A',[-s 0 0;c+s -s 0;0 d+s -s],'b',[0 1-b b],'c',[-s c d]), M0, pde, solver);
        case 'SDIRK(3,4,7)'
          s = -1.2805797612753055;
          b = 0.4453994092277531;
          c = 0.3489302860638736;
          d = 0.7586985719573739;
          e = 0.1778747841442887;
          R = RungeKuttaMethod(struct('A',[-s 0 0 0;c+s -s 0 0;0 d+s -s 0; 0 0 e+s -s],'b',[0 0 1-b b],'c',[-s c d e]), M0, pde, solver);
        % multistep BDF
        case 'BDF1'
          R = MultiStepMethod(struct('alpha',[1 -1], 'beta',[1 0]), M0, pde, solver);
        case 'BDF2'
          R = MultiStepMethod(struct('alpha',[3 -4 1], 'beta',[2 0 0]), M0, pde, solver);
        case 'BDF3'
          R = MultiStepMethod(struct('alpha',[11 -18 9 -2], 'beta',[6 0 0 0]), M0, pde, solver);
        case 'BDF4'
          R = MultiStepMethod(struct('alpha',[25 -48 36 -16 3], 'beta',[12 0 0 0 0]), M0, pde, solver);
        case 'BDF5'
          R = MultiStepMethod(struct('alpha',[137 -300 300 -200 75 -12], 'beta',[60 0 0 0 0 0]), M0, pde, solver);
        case 'BDF6'
          R = MultiStepMethod(struct('alpha',[147 -360 450 -400 225 -72 10], 'beta',[60 0 0 0 0 0 0]), M0, pde, solver);
        % multistep Adans-Bashforth
        case 'AB1'
          R = MultiStepMethod(struct('alpha',[1 -1], 'beta',[0 1]), M0, pde, solver);
        case 'AB2'
          R = MultiStepMethod(struct('alpha',[1 -1 0], 'beta',[0 3 -1]/2), M0, pde, solver);
        case 'AB3'
          R = MultiStepMethod(struct('alpha',[1 -1 0 0], 'beta',[0 23 -16 5]/12), M0, pde, solver);
        case 'AB4'
          R = MultiStepMethod(struct('alpha',[1 -1 0 0 0], 'beta',[0 55 -59 37 -9]/24), M0, pde, solver);
        case 'AB5'
          R = MultiStepMethod(struct('alpha',[1 -1 0 0 0 0], 'beta',[0 1901 -2774 2616 -1274 251]/720), M0, pde, solver);
        % multistep Adans-Moulton
        case 'AM1'
          R = MultiStepMethod(struct('alpha',[1 -1], 'beta',[1 0]), M0, pde, solver);
        case 'AM2'
          R = MultiStepMethod(struct('alpha',[1 -1], 'beta',[1 1]/2), M0, pde, solver);
        case 'AM3'
          R = MultiStepMethod(struct('alpha',[1 -1 0], 'beta',[5 8 -1]/12), M0, pde, solver);
        case 'AM4'
          R = MultiStepMethod(struct('alpha',[1 -1 0 0], 'beta',[9 19 -5 1]/24), M0, pde, solver);
        case 'AM5'
          R = MultiStepMethod(struct('alpha',[1 -1 0 0 0], 'beta',[251 646 -264 106 -19]/720), M0, pde, solver);
        case 'DG1'
          R = VariationalIntegrator(M0, pde, intType, solver);
        case 'CG2'
          R = VariationalIntegrator(M0, pde, intType, solver);
      end
    end
    function list()
      fprintf('explicit Runge Kutta Methods:\n');
      fprintf('-----------------------------\n');
      fprintf('EulerEx\n');
      fprintf('Heun2\n');
      fprintf('Heun3\n');
      fprintf('Ralston\n');
      fprintf('RK3\n');
      fprintf('RK4\n');
      fprintf('RKF\n');
      fprintf('3/8\n');
      fprintf('\n');
      fprintf('implicit Runge Kutta Methods:\n');
      fprintf('-----------------------------\n');
      fprintf('EulerIm\n');
      fprintf('CN\n');
      fprintf('Gauss2\n');
      fprintf('Gauss3\n');
      fprintf('SDIRK3\n');
      fprintf('SDIRK(3,4,5)\n');
      fprintf('SDIRK(3,4,7)\n');
      fprintf('\n');
      fprintf('Mutli Step Methods:\n');
      fprintf('-------------------\n');
      fprintf('BDF1, ... BDF6\n');
      fprintf('AB1, ... AB5\n');
      fprintf('AM1, ... AM5\n');
      fprintf('Variational Methods\n');
      fprintf('-------------------\n');
      fprintf('DG1\n');
      fprintf('CG2\n');
    end
  end
end