% PARAMETERS
dim = 2; N = 30; orderU = 1; orderV = 1; isTri = 0;
% DATA
data.a = @(x) 1 + 0*x(:,1);
data.f = @(x)sin(16*pi*prod(x,2));
dLoc = @(x)x(:,1)<Inf; 
ud = @(x)0*x(:,1);
% MESH
m = RegularMesh(N*ones(dim,1), repmat([0 1],dim,1), isTri);
if dim == 3, m.nBlock = 20; end
% FESPACES
if isTri
  fesH1 = FESpace(m, PpH1(dim,orderU), dLoc, ud);
  fesHDiv = FESpace(m, RTPp(dim,orderV));
else
  fesH1 = FESpace(m, QpH1(dim, orderU), dLoc, ud);
  fesHDiv = FESpace(m, RTQp(dim, orderV));
end
% PDE
p = DivGrad(data, fesH1, fesHDiv);
% COMPUTE
p.compute();
U = p.getSolution(1); V = p.getSolution(2);
% VISUALIZE
v = Visualizer.create(fesH1);
vHDiv = Visualizer.create(fesHDiv);
switch dim
  case 1
    v.show(U);
  case 2
    v.show(U, 'g');
  case 3
    v.show(U, 'p', struct('loc', @(x)x(:,3)<0.5));
end