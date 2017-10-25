% PARAMETERS
dim = 2; N = 100; order = 1; isTri = 1;
% DATA
clear data;
data.a = 1;
data.f = @(x)sin(16*pi*prod(x,2));
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
% FESPACE
if isTri
  e = PpL(dim, order);
else
  e = QpL(dim, order);
end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
p = Poisson(data, fes);
p.solver = IterativeSolver([], 'bicgstab', 'ilu');
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes);
switch dim
  case 1
    v.show(p.solution, 'g');
  case 2
    v.show(p.solution, 'g');
  case 3
    opt = struct('map', @(u,v)[0.5+0.5*sin(pi*v).*sin(2*pi*u), ...
                               0.5+0.5*sin(pi*v).*cos(2*pi*u), ...
                               0.5+0.5*cos(pi*v)]);
    clf; v.show(p.solution, 'g', opt);
end