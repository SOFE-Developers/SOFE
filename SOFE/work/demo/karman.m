T = 8; M = 50;
uMax = 1.5;
% DATA
clear data;
data.nu = @(x) 1e-3 + 0*x(:,1);
% MESH
m = SalomeMesh('karman');
m.uniformRefine(0);
m.nBlock = 1;
% ELEMENT
e = PpL(2,2); eP = PpL(2,1);
% FESPACE2
ud = @(x)[(4*uMax*x(:,2).*(0.41-x(:,2))/0.41^2).*(x(:,1)<eps), 0.0*x(:,1)];
fesV = FESpace(m, TPElem(e), @(x)x(:,1)<2.2, ud); 
fesP = FESpace(m, eP);
% PDE
p = EulerImplicit(Mass({1 []; [] 0}, {fesV, fesP}), ...
                  NavierStokes(data, fesV, fesP), ...
                  RegularMesh(T*M, [0, T], 0), {ud, @(x)0*x(:,1)});
% SOLVE
p.solver = DirectSolver([]);
p.integrate();
% VISUALIZE
vV = Visualizer.create(fesV);
vP = Visualizer.create(fesP);
nDoFV = fesV.getNDoF();
nDoFP = fesP.getNDoF();
opt.N = 100*[4 1]; opt.n = 20*[4 1]; opt.scale = 1; opt.width = 0.75; opt.normalize = 0;
for k = 1:5:T*M+1
  clf, vV.show(p.solution{k}(1:nDoFV), 'g', opt);
  fprintf('timestep: %d / %d\n', k, T*M+1);
  pause(0.01);
end