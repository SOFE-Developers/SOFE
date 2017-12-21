% PARAMETERS
N = 50; M = 50;
% MESH
m = RegularMesh([N; N], [0 1; 0 1], 0);
% FESPACE
fes = FESpace(m, QpL(2,1));
% PDE
p = CDR(struct('a',0.001,'b',0.5), fes);
m0 = Mass(fes);
timeline = RegularMesh(M, [0 1], 0);
u0 = @(x)0.2*exp(-((x(:,1)-0.3).^2+(x(:,2)-0.3).^2)/0.01);
% ALGORITHM
q = EulerImplicit(m0, p, timeline, u0);
q.directSolve = 2;
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
for k = 1:q.nT
  clf
  v.show(q.history{k}, 'p');
  view(3), axis([0 1 0 1 0.0 0.2]); caxis([0.0, 0.05]);
  fprintf('timestep: %d / %d\n', k, N);
  drawnow
end