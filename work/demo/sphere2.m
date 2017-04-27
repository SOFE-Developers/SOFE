% DATA
clear data;
data.a = @(x)1e4+0*x(:,1);
data.b = @(x)[0*x(:,1), 0*x(:,1), 1e3+0*x(:,1)];
% MESH
m = CADMesh([SOFE.getSOFEPath() '/meshes/library/sphere.dat'], 2);
% FESPACE
e = PpL(2, 2);
fes = FESpace(m, e);
% PDE
u0 = {@(x)exp(-(x(:,1).^2 + x(:,2).^2 + (x(:,3)-100).^2)/1000) + ...
          exp(-(x(:,1).^2 + (x(:,2)+100).^2 + x(:,3).^2)/1000)};
meshT = RegularMesh(50, [0,0.25], 0);
p = EulerImplicit(Mass(1, fes), CDR(data, fes), meshT, u0);
% SOLVE
p.integrate();
% VISUALIZE
v = Visualizer.create(fes);
N = numel(p.solution);
for k = 1:1:N
  clf, v.show(p.solution{k}, 'p', struct('n', 3)); view(-40,20)
  fprintf('timestep: %d / %d\n', k, N);
  pause(0.05);
end