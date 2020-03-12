%% PARAMETERS
N = 50; M = 100;
%% MESH
m = RegularMesh([N; N], [0 1; 0 1], 0);
%% FESPACE
fes = FESpace(m, QpL(2,1));
%% PDE
p = CDR(struct('a',0.00,'b',@(x)5*[-x(:,2)+0.5, x(:,1)-0.5]), fes);
m0 = Mass(fes);
timeline = RegularMesh(M, [0 1], 0);
u0 = @(x)0.2*exp(-((x(:,1)-0.3).^2+(x(:,2)-0.3).^2)/0.01);
%% ALGORITHM
q = TimeStep.create('dG2', Mass(fes), p, DirectSol(1));
q = Integrator(timeline, q, u0);
q.compute();
%% VISUALIZE
v = Visualizer.create(fes);
for k = 1:2:q.nT
  clf
  v.show(q.history{k}, 'p');
  view(3), axis([0 1 0 1 -0.1 0.2]); caxis([0.0, 0.1]);
  fprintf('timestep: %d / %d\n', k, M);
  drawnow
end