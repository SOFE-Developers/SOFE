%% PARAMETERS
R = 1; r = 0.5;
N = 150; isTri = 0;
M = 100; T = 2;
epsilon = 3e-2;
%% MESH
m = TorusMesh([N ceil(R/r*N)], R, r, isTri);
%% FESPACE
if isTri, e = PpL(2, 1); else, e = QpL(2, 1); end
fes = FESpace(m, e);
%% PDE
p = Poisson(struct('a',epsilon,'f',@(x,t,u)-1/epsilon*u{1}.*(u{1}.^2-1)), fes);
%% ALGORITHM
q = TimeStep.create('EulerIm', Mass(fes), p, IterativeSol('bicgstab','ilu'));
q = Integrator(RegularMesh(M, [0 T], 0), q, @(x)2*(rand(size(x,1),1)-0.5));
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
       r*sin(X), reshape(Z,size(X)));
  shading interp, axis equal
  title(sprintf('time %6.2f / %6.2f',(k-1)*T/M, T));
  pause(0.01)
end