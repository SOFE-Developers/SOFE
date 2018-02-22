clear classes
% PARAMETERS
dim = 2; N = 50; order = 1; isTri = 1;
nEig = 20;
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
m.adaptiveRefine(@(x)sum((x-[0.5 0.5]).^2,2).^0.5<0.1, 3);
m.adaptiveRefine(@(x)sum((x-[0.5 0.5]).^2,2).^0.5<0.05, 3);
m.adaptiveRefine(@(x)sum((x-[0.5 0.5]).^2,2).^0.5<0.025, 3);
%m.adaptiveRefine(@(x)sum((x-[0.5 0.5]).^2,2).^0.5<0.01, 3);
%m.show()
%return
% FESPACE
if isTri, e = PpL(dim, order); else, e = QpL(dim, order); end
fes = FESpace(m, e, @(x) x(:,1) < Inf);
% PDE
%data = struct('a',1e-6,'c',@(x)sum((x-[0.5 0.5]).^2,2));
data = struct('a',2e-1,'c',@(x)-100./sum((x-[0.5 0.5]).^2,2).^0.5);
p = CDR(data, fes);
% ALGORITHM
q = EigenSolver(p, nEig);
q.compute();
% VISUALIZE
v = Visualizer.create(fes);
disp(q.eigenVal);
for k = nEig:-1:1
  U = q.eigenVec(:,k);
  clf, v.show(U, 'l');
  axis([0 1 0 1 min(U)-0.1 max(U)+0.1]);
  pause();
end