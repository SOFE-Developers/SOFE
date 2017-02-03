% DATA
data.a    = @(x) 1+0*x(:,1);
data.f    = @(x) 1+0*x(:,1);
%% MESH
load('./meshes/library/nodesBall.dat');
load('./meshes/library/elemBall.dat');
m = Mesh(nodesBall, elemBall);
m.nBlock = 20;
% FESPACE
fes = FESpace(m, PpL(3,3), @(x)x(:,1)<Inf, @(x)0*x(:,1));
%% PDE
p = Poisson(data, fes);
p.solver = IterativeSolver([], 'bicgstab', 'ilu');
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
opt = struct('map', @(u,v)0.5*[sin(pi*v).*sin(2*pi*u), ...
                               -0.5+sin(pi*v).*cos(2*pi*u), ...
                               cos(pi*v)]);
clf; v.show(p.solution, 'g', opt);
v.show(p.solution, 'p', struct('loc', @(x)x(:,2)>0));
