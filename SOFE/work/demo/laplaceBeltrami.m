% PARAMETERS
N = 50; order = 1; isTri = 1;
% DATA
data = struct('a',0.0025,'b',[1 0 0],'f',1);
% MESH
[nodes, elem] = Mesh.getTensorProductMesh({linspace(0,1,N),linspace(0,1,N)}, isTri);
m = Mesh([nodes 0.1*sin(2*pi*nodes(:,1)).*cos(2*pi*nodes(:,2))], elem, 2);
% FESPACE
if isTri, e = PpL(2, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < eps);
% PDE
p = CDR(data, fes);
% ALGORITHM
q = IterativeSolver(p, 'bicgstab', 'ilu');
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
v.show(q.solution, 'l');