% PARAMETERS
N = 60; n = 10; K = 3;
order = 2; isTri = 0; nNonLin = 10; 
% DATA
clear data;
data.nu = 1e-3;
data.ud = @(x)[1.0*(x(:,2)>1-eps), 0.0*x(:,1)];
data.dLoc = @(x)[x(:,1)<Inf,x(:,1)<Inf];
% MESH
g1 = linspace(0,K*sqrt(data.nu),n);
g2 = linspace(K*sqrt(data.nu),1-K*sqrt(data.nu), N-2*n);
g3 = linspace(1-K*sqrt(data.nu),1,n);
grid = [g1(1:end-1) g2 g3(2:end)];
[nodes, elem] = Mesh.getTensorProductMesh({grid, grid}, isTri);
m = Mesh(nodes, elem);
% FESPACE
if isTri
  E = TPElem(PpL(2,order)); e = PpL(2,order-1);
else
  E = TPElem(QpL(2,order)); e = QpL(2,order-1);
end
FES = FESpace(m, E, data.dLoc, data.ud); fes = FESpace(m, e);
% PDE
p = NavierStokes(data, FES, fes);
% ALGORITHM
q = FixPointSolver(p, nNonLin);
q.compute();
% VISUALIZE
V = Visualizer.create(FES);
opt.scale = 0.5; opt.N = 200; opt.n = 40; opt.width = 1.5;
figure(1), V.show(p.getSolution(1), 'g', opt);
v = Visualizer.create(fes);
figure(2), v.show(p.getSolution(2), 'g', opt);
view(3), axis([0 1 0 1 -0.5 0.2]); caxis([-0.5 0.2]);