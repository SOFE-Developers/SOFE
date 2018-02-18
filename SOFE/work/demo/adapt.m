% PARAMETERS
N = 1; order = 1; isTri = 1;
% MESH
m = RegularMesh([N N], [0 1;0 1], isTri);
m.adaptiveRefine(@(x)(x(:,1)>0.2 & x(:,2)>0.2) | sum(x.^2,2).^0.5<0.2, 15);
m.coarsen(@(x)sum((x-0.5).^2,2).^0.5<0.1,15);
% FESPACE
if isTri, e = PpL(2, order); else, e = QpL(2, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
data = struct('a',1,'f', @(x)sin(50*pi*prod(x,2)));
p = Poisson(data, fes);
% ALGORITHM
q = DirectSolver(p);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
v.show(q.solution, 'g');