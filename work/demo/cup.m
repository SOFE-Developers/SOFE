% DATA
clear data;
data.nu = 1/3;
data.E = 1e6;
data.g = @(x) bsxfun(@times, [0*x(:,1) 0*x(:,1) -1e4+0*x(:,1)], (x(:,2)>80-1e-6 & x(:,3)>25));
data.dLoc = @(x) x(:,3)<-25;
% MESH;
m = CADMesh([SOFEClass.getSOFEPath() '/meshes/library/cup.dat']);
m.nBlock = 50;
% ELEMENT
e = PpL(3, 1);
% FESPACE
fes = FESpace(m, TPElem(e), data.dLoc);
% PDE
p = LinElast(data, fes);
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.show(p.solution, 'p', struct('deform', true, 'n', 1)); view(30,40);