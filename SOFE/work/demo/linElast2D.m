% PARAMETERS
N = 3; M = [5 1]'; order = 3; isTri = 0;
% DATA
data = struct('nu', 0.33,'E', 1e3, 'f', [0 -1]);
% MESH
m = RegularMesh(M*(N+1), [zeros(2,1) M], isTri);
% FESPACE
if isTri, e = PpL(2, order); else e = QpL(2, order); end
fes = FESpace(m, TPElem(e), @(x)[x(:,1)==0, x(:,1)<1]);
% PDE
p = LinElast(data, fes);
% ALGORITHM
q = DirectSolver(p);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
opt.N = M*20; opt.deform = true;
v.show(p.solution, 'g', opt);
