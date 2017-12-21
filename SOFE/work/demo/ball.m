% MESH
load([SOFE.getCorePath '/meshes/library/nodesBall.dat']);
load([SOFE.getCorePath '/meshes/library/elemBall.dat']);
% FESPACE
fes = FESpace(Mesh(nodesBall, elemBall), PpL(3,3), @(x)x(:,1)<Inf);
%% PDE
p = Poisson(struct('a',1,'f',1), fes);
% ALGORITHM
q = IterativeSolver(p, 'bicgstab', 'ilu');
q.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
opt = struct('map', @(u,v)0.5*[sin(pi*v).*sin(2*pi*u), ...
                               -0.5+sin(pi*v).*cos(2*pi*u), ...
                               cos(pi*v)]);
v.show(p.solution, 'g', opt);
v.show(p.solution, 'p', struct('loc', @(x)x(:,2)>0));
