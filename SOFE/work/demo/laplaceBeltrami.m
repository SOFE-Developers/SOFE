% PARAMETERS
N = 100; order = 1; isTri = 1;
% MESH
[nodes, elem] = Mesh.getTensorProductMesh(repmat({linspace(0,1,N)},2,1), isTri);
H = @(x) 0.5*sin(2*pi*x(:,1)).*cos(2*pi*x(:,2));
m = Mesh([nodes H(nodes)], elem, 2);
% FESPACE
if isTri, e = PpL(2, order); else, e = QpL(2, order); end
fes = FESpace(m, e, @(x) x(:,1) < 1+eps, @(x)0*x(:,1));
% PDE
data = struct('a',0.005, 'b',[1 1 0], 'f',1);
p = CDR(data, fes);
q = StaticAlg(p, DirectSol());
% SOLVE
q.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.show(q.solution);