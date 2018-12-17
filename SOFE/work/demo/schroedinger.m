% PARAMETERS
N = 50; M = 100; T = 5;
% MESH
m = RegularMesh([N; N], [0 1; 0 1], 0);
% FESPACE
fes = FESpace(m, QpL(2,1), @(x)x(:,1)<-Inf);
% PDE
p = Poisson(struct('a',0.01*i,'f',0), fes);
timeline = RegularMesh(M, [0 T], 0);
u0 = @(x)0.1+0.2*1i*exp(-((x(:,1)-0.3).^2+(x(:,2)-0.3).^2)/0.01);
% ALGORITHM
q = TimeStep.create('theta0.5', Mass(fes), p, DirectSol(1));
q = Integrator(timeline, q, u0);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
for k = 1:2:q.nT
  clf
  v.show(abs(q.history{k}).^2, 'g', struct('N',100));
  caxis([0.0, 0.03]); axis([0 1 0 1 0.0 0.08]); axis normal, view(3),
  fprintf('timestep: %d / %d\n', k-1, M);
  pause(0.05);
end