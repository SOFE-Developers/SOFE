% PARAMETERS
dim = 2; N = 1; order = 1; isTri = 1;
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
m.adaptiveRefine(@(x)(x(:,1)>0.2 & x(:,2)>0.2) | sum(x.^2,2).^0.5<0.2, 15);
m.coarsen(@(x)sum((x-0.5).^2,2).^0.5<0.1,15);
% FESPACE
if isTri, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
data = struct('a',1,'f', @(x)sin(50*pi*prod(x,2)));
p = Poisson(data, fes);
p.createSys = 1;
% ALGORITHM
q = IterativeSolver(p, 'bicgstab', 'ilu');
%q = DirectSolver(p);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
switch dim
  case 1
    v.show(q.solution, 'g');
  case 2
    v.show(q.solution, 'g');
  case 3
    opt = struct('map', @(u,v)[0.5+0.5*sin(pi*v).*sin(2*pi*u), ...
                               0.5+0.5*sin(pi*v).*cos(2*pi*u), ...
                               0.5+0.5*cos(pi*v)]);
    v.show(q.solution, 'g', opt);
end