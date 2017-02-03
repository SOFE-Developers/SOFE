% PARAMETERS
N = 10; M = [5 1 1]; order = 1; isTri = 1;
% DATA
data.nu = 1/3; data.E = 1e3;
data.f = @(x) bsxfun(@plus, [0 0 -1], 0*x(:,1));
% MESH;
m = RegularMesh(M'*(N+1), [zeros(1,3);M]', isTri);
m.nBlock = min(20, m.topology.getNumber(3));
% ELEMENT
if isTri, e = PpL(3, order); else  e = QpS(3, order); end
% FESPACE
fes = FESpace(m, TPElem(e), @(x)x(:,1)==0);
% PDE
p = LinElast(data, fes);
p.solver = IterativeSolver([], 'bicgstab', 'ilu');% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
opt.N = 50;;
opt.scale = 1.0; opt.n = 50; opt.width = 2.0;
opt.deform = 1;
opt.map = @(u,v)[M(1)*u, 0.5+0.5*M(2)*cos(2*pi*v), 0.5+0.5*M(3)*sin(2*pi*v)];
v.show(p.solution, 'g', opt);
