% PARAMETERS
order = 12;
% MESH
load([SOFE.getCorePath '/meshes/library/nodesRing.dat']);
load([SOFE.getCorePath '/meshes/library/elemRing.dat']);
m = Mesh(nodesRing, elemRing);
% FESPACE
fes = FESpace(m, PpL(2,order), @(x)sum(x.^2,2).^0.5<0.5 & x(:,1)<0);
% PDE
p = Poisson(struct('a',1,'f',1), fes);
% ALGORITHM
q = IterativeSolver(p, 'bicgstab', 'ilu');
q.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.show(q.solution, '', struct('n',order));
