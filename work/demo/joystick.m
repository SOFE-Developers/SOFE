% DATA
data.nu = 1/3;
data.E = 1e6;
data.f = @(x) bsxfun(@plus, [0 0 -0.0], 0*x(:,1));
data.g = @(x) bsxfun(@times, 2*[200+0*x(:,1) 0*x(:,1) 0*x(:,1)], (x(:,3)>450));
data.dLoc = @(x)x(:,3) == 0;
% MESH;
m = CADMesh([SOFEClass.getSOFEPath() '/meshes/library/joystick.dat']);
% ELEMENT
e = PpL(3, 2);
% FESPACE
fes = FESpace(m, TPElem(e), data.dLoc);
% PDE
p = LinElast(data, fes);
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.show(p.solution, 'p', struct('deform', true, 'n', 1)); view(30,40);
%
%hold on, m.show, hold off