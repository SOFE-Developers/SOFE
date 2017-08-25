N = 60; nNonLin = 7;
isTri = 1; order = 2;
% DATA
clear data;
data.nu = @(x) 1e-3 + 0*x(:,1);
data.ud = @(x)[1.0*(x(:,2)>1-eps), 0.0*x(:,1)];
data.dLoc = @(x)[x(:,1)<Inf,x(:,1)<Inf];
% MESH
m = RegularMesh([N N], [0 1; 0 1], isTri);
% ELEMENT
if isTri
  e = PpL(2,order); eP = PpL(2,order-1);
else
  e = QpL(2,order); eP = QpL(2,order-1);
end
% FESPACE
fesV = FESpace(m, TPElem(e), data.dLoc, data.ud);
fesP = FESpace(m, eP);
% PDE
p = NavierStokes(data, fesV, fesP);
% SOLVE % VISUALIZE
vV = Visualizer.create(fesV);
vP = Visualizer.create(fesP);
opt.scale = 0.5; opt.N = 200; opt.n = 40; opt.width = 1.5;
for k = 1:nNonLin
  p.compute();
  p.setState(p.solution);
  figure(1), vV.show(p.getSolution(1), 'g', opt);
  fprintf('step: %d / %d\n', k, nNonLin);
  pause(0.01);
end