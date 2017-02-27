% PARAMETERS
dim = 1; isTri = 0;
T = 10; N = 80; M = 50;
center = 0.5*ones(1,dim);
% DATA
data.a = @(x) 0.05 + 0*x(:,1);
data.f = @(x) 0*x(:,1);
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), 0);
% FESPACES
fes = FESpace(m, QpL(dim,1));
% TIME
u0 = @(x)0.1*exp(-sum(bsxfun(@minus, x, center).^2,2)/0.005);
meshT = RegularMesh(T*M, [0 1], 0); meshT.topology.scale(T);
% PDE
pTime = LeapFrog(Mass(1, fes), Poisson(data, fes), meshT, u0);
% COMPUTE
pTime.integrate();
% VISUALIZE
v = Visualizer.create(fes);
N = numel(pTime.solution);
for k = 1:5:N
  clf
  switch dim
    case 1
      v.show(pTime.solution{k}, 0);
      axis([0 1.0 -0.1 0.1]);
    case 2
      v.show(pTime.solution{k}, 'g');
      caxis([-0.01, 0.01]);
  end
  pause(0.05);
  fprintf('timestep: %d / %d\n', k, N);
end