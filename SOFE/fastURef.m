% PARAMETERS
dim = 2; N = 40; order = 1; isTri = 1;
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
%tic, m.uniformRefine(2); toc

%return
%m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
%tic,
%m.uniformRefineFast(2);
%toc
%m.topology.getNumber(2)

% FESPACE
if isTri, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
p = Poisson(struct('a',1,'f', @(x)sin(16*pi*prod(x,2))), fes);
p.createSys = 1;
% ALGORITHM
q = IterativeSolver(p, 'bicgstab', 'none');
q.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
if dim==2
  v.show(q.solution, 'p', struct('n', order));
else
  v.show(q.solution, 'g',struct('map',@(u,v)[u v 0.5+0*u]));
end