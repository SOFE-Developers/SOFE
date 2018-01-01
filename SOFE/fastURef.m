% PARAMETERS
dim = 2; N = 100; order = 1; isTri = 1;
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
tic,
m.uniformRefine(1);
toc
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
tic,
m.uniformRefineFast(1);
toc
% FESPACE
if isTri, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
p = Poisson(struct('a',1,'f', @(x)sin(16*pi*prod(x,2))), fes);
p.createSys = 1;
% ALGORITHM
q = IterativeSolver(p, 'bicgstab', 'ilu');
q.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.show(q.solution, 'p');