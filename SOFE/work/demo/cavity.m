% PARAMETERS
N = 50; n = 10; K = 5;
order = 2; isTri = 0; nNonLin = 10; 
% DATA
data = struct('nu', 1e-4, 'f', 0);
data.ud = @(x)[(x(:,1).*(1-x(:,1))).*(x(:,2)>1-eps), 0.0*x(:,1)];
%data.ud = @(x)[x(:,1).*(1-x(:,1)).*(x(:,1)-0.5).*(x(:,2)>1-eps), 0.0*x(:,1)];
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
FES = FESpace(m, E, data.dLoc, data.ud);
fes = FESpace(m, e);
% PDE
p = NavierStokes(data, FES, fes);
% ALGORITHM
q = FixPointSolver(p, nNonLin);
q.compute();
% VISUALIZE
V = Visualizer.create(FES);
opt.scale = 0.5; opt.N = 200; opt.n = 40; opt.width = 1.5;
figure(1), clf
U = q.linSolver.solution(p.J(1,1):p.J(1,2));
V.show(U, 'g', opt);
v = Visualizer.create(fes);
P = q.linSolver.solution(p.J(2,1):p.J(2,2));
figure(2), clf
v.show(P, 'g', opt); view(3),axis normal