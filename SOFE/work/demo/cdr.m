% PARAMETERS
dim = 2; N = 30; order = 1; isTri = 0;
% DATA
clear data;
data.a = @(x)0.1+0*x(:,1);
data.b = @(x)1+0*x;
data.c = @(x)1+0*x(:,1);
data.f = @(x)1+0*x(:,1);
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
if dim == 3, m.nBlock = 20; end
% FESPACE
if isTri
  e = PpL(dim, order);
else
  e = QpL(dim, order);
end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
p = CDR(data, fes);
if dim == 3
  p.solver = IterativeSolver([], 'bicgstab', 'ilu');
end
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
switch dim
  case 1
    v.show(p.solution);
  case 2
    v.show(p.solution, 'g');
  case 3
    opt = struct('map', @(u,v)[0.5+0.5*sin(pi*v).*sin(2*pi*u), ...
                               0.5+0.5*sin(pi*v).*cos(2*pi*u), ...
                               0.5+0.5*cos(pi*v)]);
    v.show(p.solution, 'g', opt);
end