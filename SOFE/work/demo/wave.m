% PARAMETERS
dim = 1; T = 10; N = 80; M = 500; isTri = 1;
cc = 0.5*ones(1,dim);
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), 0);
% FESPACE
fes = FESpace(m, QpL(dim,1),@(x)x(:,1)<Inf);
% PDE
p = Poisson(struct('a',0.05,'f',0), fes);
m0 = Mass(fes);
timeline = RegularMesh(M, [0 1], 0); timeline.scale(T);
u0 = @(x)0.1*exp(-sum(bsxfun(@minus, x, cc).^2,2)/0.005);
% ALGORITHM
q = LeapFrog(m0, p, timeline, u0);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
for k = 1:5:q.nT
  clf
  switch dim
    case 1
      v.show(q.history{k}, 0); axis([0 1.0 -0.1 0.1]);
    case 2
      v.show(q.history{k}, 'g');caxis([-0.01, 0.01]);
  end
  drawnow
  fprintf('timestep: %d / %d\n', k, q.nT);
end