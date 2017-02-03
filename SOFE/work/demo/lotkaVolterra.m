N = 30; T = 5; M = 50;
% DATA
clear data;
alpha = 1; beta = 8/3; gamma = 1; delta = 1;
data.a = {@(x) 0.01 + 0*x(:,1), @(x) 0.01 + 0*x(:,1)};
data.f{1} = @(x, t, U) U{1}.*(alpha - beta*U{2});
data.f{2} = @(x, t, U) -U{2}.*(gamma - delta*U{1});
% MESH
m = RegularMesh([N; N], [0 1;0 1], 0);
% FESPACES
fes = FESpace(m, QpH1(2,1), @(x)x(:,1)<-Inf);
% time
meshT = RegularMesh(T*M, [0 T], 0);
u0 = {@(x)5*exp(-((x(:,1)-0.4).^2+(x(:,2)-0.4).^2)/0.01), ...
      @(x)5*exp(-((x(:,1)-0.6).^2+(x(:,2)-0.6).^2)/0.01)};
% PROBLEM
p = EulerImplicit(Mass(eye(2), {fes, fes}), LotkaVolterra(data, fes), meshT, u0);
p.solver = DirectSolver([], 1);
% COMPUTE
p.integrate();
% VISUALIZE
v = Visualizer.create(fes);
N = numel(p.solution);
for k = 1:5:N
  sol1 = p.solution{k}(1:end/2);
  sol2 = p.solution{k}(end/2+1:end);
  figure(1), clf, v.patch(sol1);
  figure(2), clf, v.patch(sol2);
  fprintf('timestep: %d / %d\n', k, N);
  drawnow;
end