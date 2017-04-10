% DATA
clear data;
data.a = @(x)1e5+0*x(:,1);
data.b = @(x)[50*x(:,2), -1e2*x(:,1), 0*x(:,1)];
% MESH
m = CADMesh([SOFEClass.getSOFEPath() '/meshes/library/torus.dat'], 2);
% FESPACE
e = PpL(2, 1);
fes = FESpace(m, e);
% PDE
u0 = {@(x)exp(-(x(:,1).^2 + (x(:,2)+400).^2 + x(:,3).^2)/10000)};
meshT = RegularMesh(50, [0,0.25], 0);
p = EulerImplicit(Mass(1, fes), CDR(data, fes), meshT, u0);
% SOLVE
p.integrate();
% VISUALIZE
v = Visualizer.create(fes);
N = numel(p.solution);
for k = 1:1:N
  clf, v.show(p.solution{k}, 'p', struct('n', 2)); view(-40,20), axis equal
  fprintf('timestep: %d / %d\n', k, N);
  pause(0.01);
end