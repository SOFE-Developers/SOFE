% PARAMETERS
dim = 2; N = 50; order = 1; isTri = 1;
% DATA
clear data;
data.a = @(x)0.0025+0*x(:,1);
data.b = @(x)1*[x(:,1) 0*x(:,1) 0*x(:,1)];
data.f = @(x)1+0*sin(4*pi*x(:,1));
% MESH
m = [];
[nodes, elem] = Mesh.getTensorProductMesh(repmat({linspace(0,1,N)},dim,1), isTri);
if dim == 1
  m = Mesh([nodes sin(2*pi*nodes) cos(2*pi*nodes)], elem, dim);
elseif dim == 2
  m = Mesh([nodes sin(2*pi*nodes(:,1)).*cos(2*pi*nodes(:,2))], elem, dim);
end
% FESPACE
if isTri
  e = PpL(dim, order);
else
  e = QpL(dim, order);
end
fes = FESpace(m, e, @(x) x(:,1) < eps, @(x)0*x(:,1));
% PDE
p = CDR(data, fes);
% SOLVE
p.compute();
% VISUALIZE
v = Visualizer.create(fes); clf
v.show(p.solution, 'l', struct('N',300));