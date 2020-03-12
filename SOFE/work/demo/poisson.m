%% PARAMETERS
dim = 2; N = 30; order = 2; isTri = 1;
%% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
%% FESPACE
if isTri==1, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
%% PDE
data = struct('a',1,'f', @(x)sin(16*pi*prod(x,2)));
p = Poisson(data, fes);
%% SOLVE
q = StaticAlg(p, IterativeSol('bicgstab', 'ilu'));
% q = StaticAlg(p, DirectSol());
q.compute();
%% VISUALIZE
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