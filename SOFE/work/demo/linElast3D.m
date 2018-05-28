% PARAMETERS
N = 10; M = [5 1 1]; order = 1; isTri = 1;
% DATA
data = struct('nu', 0.33,'E', 1e3, 'f', [0 0 -1]);
% MESH
m = RegularMesh(M'*(N+1), [zeros(1,3);M]', isTri);
% FESPACE
if isTri, e = TPElem(PpL(3, order)); else  e = TPElem(QpL(3, order)); end
fes = FESpace(m, e, @(x)[x(:,1)==0, x(:,1)<1, x(:,1)<1]);
% PDE
p = LinElast(data, fes);
% ALGORITHM
q = StaticAlg(p, IterativeSol('bicgstab','ilu'));
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
opt.N = 50; opt.deform = true;
opt.map = @(u,v)[M(1)*u, 0.5+0.5*M(2)*cos(2*pi*v), 0.5+0.5*M(3)*sin(2*pi*v)];
v.show(q.solution, 'g', opt);
