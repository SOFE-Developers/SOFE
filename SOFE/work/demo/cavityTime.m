% PARAMETERS
N = 50; n = 10; K = 4;
M = 50; T = 3e2;
order = 2; isTri = 0;
% DATA
data = struct('nu', 1e-5, 'f', 0);
data.ud = @(x)[(x(:,1).*(1-x(:,1))).*(x(:,2)>1-eps), 0.0*x(:,1)];
data.dLoc = @(x)[x(:,1)<Inf,x(:,1)<Inf];
% MESH
n = 10; K = 5;
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
%M0 = Mass2({1 {} {}; {} 0 {}; {} {} 0},{FES, fes ScalarFESpace()});
M0 = Mass2({1 {}; {} 0},{FES, fes});
timeline = RegularMesh(M, [0 T], 0);
u0 = @(x)0.2*exp(-((x(:,1)-0.3).^2+(x(:,2)-0.3).^2)/0.01);
q = EulerImplicit(M0, p, timeline); q.directSolve = 1;
q.compute();
% VISUALIZE
V = Visualizer.create(FES);
opt = struct('abs',1,'vectors',1,'normalize',0,'scale',1.0);
for k = 1:q.nT
  clf
  V.show(q.history{k}(p.J(1,1):p.J(1,2)), 'g', opt);
  fprintf('timestep: %d / %d\n', k, N);
  drawnow
end