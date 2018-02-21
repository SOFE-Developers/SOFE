clear classes
% PARAMETERS
dim = 2; N = 20; order = 1; isTri = 2;
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
%m.adaptiveRefine(@(x)x(:,1)<0.75 & x(:,2)<0.5,5);
% FESPACE
if isTri, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
data = struct('a',1,'c',-0);
p = CDR(data, fes);
% ALGORITHM
q = EigenSolver(p);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
disp(q.eigenVal);
v.show(q.eigenVec(:,16), 'g');
axis tight