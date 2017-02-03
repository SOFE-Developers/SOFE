N = 50; nNonLin = 7;
isTri = 1; order = 2;
% DATA
data.nu = @(x) 1e-3 + 0*x(:,1);
data.ud = @(x)[1.0*(x(:,2)>1-eps), 0.0*x(:,1)];
data.dLoc = @(x)[x(:,1)<Inf,x(:,1)<Inf];
% MESH
m = RegularMesh([N N], [0 1; 0 1], isTri);
m.nBlock = 1;
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
opt.scale = 0.5; opt.N = 200; opt.n = 20; opt.width = 1.5;
for k = 1:nNonLin
  p.compute();
  figure(1), vV.show(p.getSolution(1), 'g', opt);
  fprintf('step: %d / %d\n', k, nNonLin);
  pause(0.01);
end
figure(2), vP.surf(p.getSolution(2)); axis([0 1 0 1]); view(-30,30);