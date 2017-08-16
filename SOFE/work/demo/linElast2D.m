% PARAMETERS
N = 3; M = [5 1]'; order = 3; isTri = 0;
% DATA
clear data;
data.nu = 1/3; data.E = 1e3;
data.f = @(x) bsxfun(@plus, [0 -1], 0*x(:,1));
% MESH
m = RegularMesh(M*(N+1), [zeros(2,1) M], isTri);
% ELEMENT
if isTri, e = PpL(2, order); else e = QpL(2, order); end
% FESPACE
%fes = FESpace(m, TPElem(e), @(x)[x(:,1)==0, x(:,1)<2], @(x)0*x);
fes = FESpace(m, TPElem(e), @(x)x(:,1)==0, @(x)0*x);
% PDE
p = LinElast(data, fes);
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
opt.N = M*20; opt.n = 50; opt.deform = true; opt.scale = 1.0; opt.width = 2;
v.show(p.solution, 'g', opt);
