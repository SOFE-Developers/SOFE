%% PARAMETERS
dim = 3; N = 3; ll = [5 ones(1,dim-1)]'; order = 1; isTri = 1;
T = 5; M = 4000;
%% MESH
m = RegularMesh(ll*(N+1), [zeros(dim,1) ll], isTri);
%% ELEMENT
if isTri, e = PpL(dim, order); else e = QpL(dim, order); end
%% FESPACE
fes = FESpace(m, VecElem(e), @(x)x(:,1)==0, @(x)0*x);
%% PDE
p = LinElast(struct('nu', 0.33,'E', 1e3, 'f', [zeros(1,dim-1) -1]), fes);
m0 = Mass(fes);
timeline = RegularMesh(M, [0 1], 0); timeline.scale(T);
u0 = @(x)0*x;
%% ALGORITHM
q = LeapFrog(m0, p, DirectSol(1));
q = Integrator(timeline, q, u0);
q.compute();
%% VISUALIZE
v = Visualizer.create(fes); clf
opt.N = ll*20; opt.deform = true;
if dim == 3
  opt.map = @(u,v)[ll(1)*u, 0.5+0.5*ll(2)*cos(2*pi*v), 0.5+0.5*ll(3)*sin(2*pi*v)];
  ax = [0 5.5 0 1 -2 1.5];
else
  ax = [0 5.5 -2 1.5];
end
for k = 1:100:q.nT
  clf, v.show(0.2*q.history{k}, 'g', opt);
  axis(ax); caxis([[0 1]]);
  drawnow
  fprintf('timestep: %d / %d\n', k, q.nT);
end
