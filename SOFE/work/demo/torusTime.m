%% PLUGINS
SOFE.plugin('integrator');
%% PARAMETERS
R = 1; r = 0.25;
N = 20; order = 3; isTri = 0;
M = 100; T = 15;
%% MESH
m = TorusMesh([N ceil(R/r*N)], R, r, isTri);
%% FESPACE
if isTri, e = PpL(2,order); else e = QpL(2,order); end
fes = FESpace(m, e);
%% PDE
data = struct('a',0.001, 'b', @(x)0.5*[-x(:,2), x(:,1), 0*x(:,1)] + ...
      [-x(:,3).*x(:,1)./sqrt(x(:,1).^2 + x(:,2).^2), -x(:,3).*x(:,2)./sqrt(x(:,1).^2 + x(:,2).^2), sqrt(x(:,1).^2 + x(:,2).^2)-R]);
initialC = @(x)0.2*exp(-(sum((x-[-(R+r) 0 0]).^2,2))/0.1);
p = CDR(data, fes);
q = TimeStep.create('DG1', Mass(fes), p, DirectSol(1));
q = Integrator(RegularMesh(M, [0 T]), q, initialC);
q.compute();
%% VISUALIZE
v = Visualizer.create(fes);
mm = RegularMesh([N ceil(R/r*N)], [0 2*pi;0 2*pi], isTri);
[X,Y] = meshgrid(linspace(0,2*pi,200));
P = [X(:) Y(:)];
P = mm.evalInversReferenceMap(P);
for k = 1:1:q.nT
  Z = fes.evalDoFVector(q.history{k}, P, [], 0);
  surf((R+r*cos(X)).*cos(Y), ...
       (R+r*cos(X)).*sin(Y), ...
       r*sin(X), ...
       reshape(Z,size(X)));
  shading interp, axis equal
  title(sprintf('time %6.2f / %6.2f',(k-1)*T/M, T));
  pause(0.02)
end
