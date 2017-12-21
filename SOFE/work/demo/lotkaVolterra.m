N = 30; T = 10; M = 20;
% DATA
clear data;
alpha = 1; beta = 8/3; gamma = 1; delta = 1;
data.a = {@(x) 0.01 + 0*x(:,1), @(x) 0.01 + 0*x(:,1)};
data.f{1} = @(x, t, U) U{1}.*(alpha - beta*U{2});
data.f{2} = @(x, t, U) -U{2}.*(gamma - delta*U{1});
% MESH
m = RegularMesh([N; N], [0 1;0 1], 0);
% FESPACE
fes = FESpace(m, QpL(2,1));
% PDE
p = LotkaVolterra(data, fes);
m0 = Mass(fes, 2);
timeline = RegularMesh(T*M, [0 T], 0);
u0 = {@(x)5*exp(-((x(:,1)-0.4).^2+(x(:,2)-0.4).^2)/0.01), ...
      @(x)5*exp(-((x(:,1)-0.6).^2+(x(:,2)-0.6).^2)/0.01)};
% ALGORITHM
q = EulerImplicit(m0, p, timeline, u0);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
for k = 1:2:q.nT
  U1 = q.history{k}(1:end/2);
  U2 = q.history{k}(end/2+1:end);
  figure(1), clf, v.show(U1, 'p'); axis([0 1 0 1 -10 10]); caxis([0 3]);
  figure(2), clf, v.show(U2, 'p'); axis([0 1 0 1 -10 10]); caxis([0 1]);
  fprintf('timestep: %d / %d\n', k, T*M);
  drawnow;
end