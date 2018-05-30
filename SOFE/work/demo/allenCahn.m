%% PARAMETERS
N = 200; M = 100; epsilon = 2e-2; T = 1;
%% MESH
m = TorusMesh([N N], 1, 0.5, 0);
%% FESPACE
fes = FESpace(m, QpL(2,1));
%% PDE
p = Poisson(struct('a',epsilon,'f',@(x,t,u)-1/epsilon*u{1}.*(u{1}.^2-1)), fes);
%% ALGORITHM
q = TimeStep.create('theta0.5', Mass(fes), p, IterativeSol('bicgstab','ilu'));
q = Integrator(RegularMesh(M, [0 T], 0), q, @(x)2*(rand(size(x,1),1)-0.5));
q.compute();
%% VISUALIZE
v = Visualizer.create(fes);
for k = 1:2:q.nT
  clf
  v.show(q.history{k}, 'p');
  view(3),
  fprintf('timestep: %d / %d\n', k, q.nT);
  title(sprintf('time %6.2f',(k-1)*T/M));
  drawnow
end