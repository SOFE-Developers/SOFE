% PARAMETERS
dim = 3; N = 3; M = [5 ones(1,dim-1)]'; order = 2; isTri = 0;
% MESH
m = RegularMesh(M*(N+1), [zeros(dim,1) M], isTri);
% FESPACE
if isTri, e = PpL(dim, order); else e = QpL(dim, order); end
fes = FESpace(m, TPElem(e), @(x)x(:,1)==0, [0 0 0]);
% PDE
p = LinElast(struct('nu', 0.33,'E', 1e3, 'f', [zeros(1,dim-1) -1]), fes);
q = StaticAlg(p, DirectSol());
q.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
opt.N = M*30; opt.deform = true;
if dim == 3
  opt.map = @(u,v)[M(1)*u, 0.5+0.5*M(2)*cos(2*pi*v), 0.5+0.5*M(3)*sin(2*pi*v)];
end
v.show(q.solution, 'g', opt);