N = 50; M = 50;
% DATA
clear data;
data.a    = @(x) 0.01 + 0*x(:,1);
data.b    = @(x) -0.5*[1+0*x(:,1) 1+0*x(:,1)];
data.g    = @(x) 0*x(:,1);
% MESH
m = RegularMesh([N; N], [0 1; 0 1], 0);
% FESPACES
fes = FESpace(m, QpL(2,1), @(x)x(:,1)>0.5, @(x)0.1*(sum(x,2)>1.8));
% TIME
u0 = @(x)0.2*exp(-((x(:,1)-0.3).^2+(x(:,2)-0.3).^2)/0.01);
meshT = RegularMesh(M, [0 1], 0);
% PDE
p = CDR2(data, fes); p.createSys = 0;
q = EulerImplicit2(Mass2(fes), p, meshT, u0);
% COMPUTE
q.integrate();
% VISUALIZE
v = Visualizer2D(fes);
for k = 1:1:q.nT
  clf
  v.show(q.history{k}, 'p');
  view(3), axis([0 1 0 1 0.0 0.2])
  caxis([0.0, 0.05]);
  fprintf('timestep: %d / %d\n', k, N);
  pause(0.05);
end