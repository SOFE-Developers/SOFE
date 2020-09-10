clear classes
%% PARAMETERS
N = 20; order = 2; isTri = 1;
%% MESH
m = RegularMesh([N N], [0 1; 0 1], isTri);
%% FESPACE
if isTri
  E = VecElem(PpL(2,order+1)); e = PpL(2,order);
else
  E = VecElem(QpL(2,order+1)); e = QpL(2,order);
end
FES = FESpace(m, E, @(x)x < Inf);
fes = FESpace(m, e);
%% PDE
data = struct('nu', 1, 'f', @(x)[-x(:,2)+0.5 x(:,1)-0.5]);
p = Stokes(data, FES, fes);
%% ALGORITHM
q = StaticAlg(p, DirectSol());
q.assemble();
fes.getFreeDoFs();
fes.freeDoFs(1) = false;
q.solve()
%% VISUALIZE
V = Visualizer.create(FES);
v = Visualizer.create(fes);
opt.scale = 0.5; opt.N = 200; opt.n = 40; opt.width = 1.5;
figure(1), clf
U = q.solution(p.J(1,1):p.J(1,2)); V.show(U, 'g', opt);
%%P = q.solution(p.J(2,1):p.J(2,2)); v.show(P, 'g', opt); view(3); axis normal
drawnow
norm(U)