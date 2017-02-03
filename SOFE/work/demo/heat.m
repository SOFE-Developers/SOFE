N = 50; M = 50;
% DATA
clear data;
data.a    = @(x) 0.01 + 0*x(:,1);
data.b    = @(x)0.5*[1+0*x(:,1) 1+0*x(:,1)];;
data.f    = @(x) 0*x(:,1);
data.g    = @(x) 0*x(:,1);
% MESH
m = RegularMesh([N; N], [0 1; 0 1], 0);
% FESPACES
fes = FESpace(m, QpH1(2,1));
% TIME
u0 = {@(x)0.2*exp(-((x(:,1)-0.3).^2+(x(:,2)-0.3).^2)/0.01)};
meshT = RegularMesh(M, [0 1], 0);
% PDE
p = EulerImplicit(Mass(1, fes), CDR(data, fes), meshT, u0);
p.solver = DirectSolver([], 1);
% COMPUTE
p.integrate();
% VISUALIZE
v = Visualizer2D(fes);
N = numel(p.solution);
for k = 1:1:N
  clf
  v.show(p.solution{k}, 'p');
  view(-60,40), axis([0 1 0 1 0.0 0.2])
  caxis([0.0, 0.05]);
  fprintf('timestep: %d / %d\n', k, N);
  pause(0.05);
end