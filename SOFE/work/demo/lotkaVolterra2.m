N = 30; T = 10; M = 25;
% DATA
clear data;
alpha = 1; beta = 8/3; gamma = 1; delta = 1;
data.a = {@(x) 0.01 + 0*x(:,1), @(x) 0.01 + 0*x(:,1)};
data.f{1} = @(x, t, U) U{1}.*(alpha - beta*U{2});
data.f{2} = @(x, t, U) -U{2}.*(gamma - delta*U{1});
% MESH
m = RegularMesh([N; N], [0 1;0 1], 0);
% FESPACES
fes = FESpace(m, QpL(2,1), @(x)x(:,1)<-Inf);
% time
meshT = RegularMesh(T*M, [0 T], 0);
u0 = {@(x)5*exp(-((x(:,1)-0.4).^2+(x(:,2)-0.4).^2)/0.01), ...
      @(x)5*exp(-((x(:,1)-0.6).^2+(x(:,2)-0.6).^2)/0.01)};
% PROBLEM
p = LotkaVolterra2(data, fes); p.createSys = 0;
q = EulerImplicit2(Mass2(fes, 2), p, meshT, u0);
% COMPUTE
q.integrate();
% VISUALIZE
v = Visualizer.create(fes);
for k = 1:5:q.nT
  sol1 = q.history{k}(1:end/2);
  sol2 = q.history{k}(end/2+1:end);
  figure(1), clf, v.patch(sol1);
  figure(2), clf, v.patch(sol2);
  fprintf('timestep: %d / %d\n', k, N);
  drawnow;
end