%% PARAMETERS
dim = 3; N = 20; order = 1; isTri = 1;
%% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
%% FESPACE
if isTri, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
%% PDE
p = Poisson(struct('a',1,'f',@(x)sin(16*pi*prod(x,2))), fes);
q = EigenSolver(p, 50);
q.compute(); 
%% VISUALIZE
v = Visualizer.create(fes); clf
disp(q.eigenVal)
for k = 1:numel(q.eigenVal)
  U = q.eigenVec(:,k);
  switch dim
    case 1
      v.show(U, 'g');
    case 2
      v.show(U, 'g');
    case 3
      v.show(U, 'g', struct('map', @(u,v)[u, v, 0.5+0*v])); hold on
      v.show(U, 'g', struct('map', @(u,v)[0.5+0*v, u, v]));
      v.show(U, 'g', struct('map', @(u,v)[u, 0.5+0*v, v])); hold off
  end
  pause()
end