N = 50; M = 30; T = 3;
% DATA
data.a    = @(x) -1i*0.05 + 0*x(:,1);
data.f    = @(x) 0*x(:,1);
data.dLoc = @(x) x(:,1) < Inf;
data.ud   = @(x) 0*x(:,1);
% MESH
m = RegularMesh([N; N], [0 1; 0 1], 0);
% FESPACES
fes = FESpace(m, QpH1(2,1), data.dLoc);
% TIME
u0 = {@(x)0.1+0.2*1i*exp(-((x(:,1)-0.3).^2+(x(:,2)-0.3).^2)/0.01)};
meshT = RegularMesh(T*M, [0 T], 0);
% PDE
p = EulerImplicit(Mass(1, fes), Poisson(data, fes), meshT, u0);
p.solver = DirectSolver([], 1);
% COMPUTE
p.integrate();
% VISUALIZE
v = Visualizer.create(fes);
nT = numel(p.solution);
for k = 1:2:nT
  clf
  v.show(abs(p.solution{k}).^2, 'g', struct('N',100));
  caxis([0.0, 0.03]); axis([0 1 0 1 0.0 0.1]); axis normal, view(-60,40),
  fprintf('timestep: %d / %d\n', k, N);
  pause(0.05);
end