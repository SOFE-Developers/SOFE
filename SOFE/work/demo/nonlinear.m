% PARAMETERS
dim = 2; N = 50; order = 1; isTri = 1;
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
% FESPACE
if isTri, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
data = struct('a',3e-2,'b',@(x,t,u,du)x(:,2).*du{1},'f',@(x)x(:,1));
p = CDR(data, fes);
% SOLVE
q = FixPoint(p, 30, DirectSol());
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
figure(1), clf
for k = 1:nNonLin
  v.show(q.history{k}, 'g', opt); view(3); axis normal
  title(sprintf('Iteration: %d / %d', k, nNonLin));
  drawnow
end