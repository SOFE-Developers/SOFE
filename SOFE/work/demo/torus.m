%% MESH
m = TorusMesh([30 40], 1, 0.2, 1);
%% FESPACE
fes = FESpace(m, PpL(2,2));
%% PDE
p = CDR(struct('a',1, 'c',@(x)x(:,1), 'f',1), fes);
%% ALGORITHM
q = StaticAlg(p, DirectSol());
q.compute();
%% VISUALIZE
v = Visualizer.create(fes); clf
v.patch(q.solution, struct('n',2));
view(3)
