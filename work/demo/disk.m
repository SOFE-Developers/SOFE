% PARAMETERS
dim = 3; order = 1; isTri = 1;
% DATA
clear data;
data.a = 1;
data.f = @(x)1+0*x(:,1);
% MESH
m = SalomeMesh('disk'); m.uniformRefine(1);
m.nBlock = 20;
% FESPACE
if isTri
  e = PpL(dim, order);
else
  e = Qp(dim, order);
end
fes = FESpace(m, e, @(x) x(:,1) < Inf, []);
% PDE
p = Poisson(data, fes);
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.show(p.solution,'g',struct('map',@(u,v)[0*u, -100+200*u, -50+70*v]));
hold on
m.show();
hold off