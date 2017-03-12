% DATA
data.a    = @(x) 1+0*x(:,1);
data.f    = @(x) 1+0*x(:,1);
%% MESH
load([SOFEClass.getSOFEPath '/meshes/library/nodesRing.dat']);
load([SOFEClass.getSOFEPath '/meshes/library/elemRing.dat']);
m = Mesh(nodesRing, elemRing);
% FESPACE
fes = FESpace(m, PpL(2,5), @(x)sum(x.^2,2).^0.5<10.5, @(x)0*x(:,1));
%% PDE
p = Poisson(data, fes);
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.patch(p.solution, struct('n',10)); axis normal
