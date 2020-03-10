%% PLUGINS
SOFE.plugin('integrator');
%% PARAMETERS
N = 50; T = 10; M = 20;
%% DATA
clear data;
alpha = 1; beta = 8/3; gamma = 1; delta = 1;
data.a = {@(x) 1e-3 + 0*x(:,1), @(x) 5e-3 + 0*x(:,1)};
data.f{1} = @(x, t, U) U{1}.*(alpha.*(1 - U{1}/1e3) - beta*U{2});
data.f{2} = @(x, t, U) -U{2}.*(gamma - delta*U{1});
%% PDE
m = RegularMesh([N; N], [0 1;0 1], 1);
fes = FESpace(m, PpL(2,1), @(x)x(:,1)<-Inf);
p = LotkaVolterra(data, fes);
timeline = RegularMesh(T*M, [0 T], 0);
%
initialC = {@(x)5*exp(-((x(:,1)-0.4).^2+(x(:,2)-0.4).^2)/0.01), ...
            @(x)5*exp(-((x(:,1)-0.6).^2+(x(:,2)-0.6).^2)/0.01)};
q = Integrator(timeline, TimeStep.create('theta0.5', Mass(fes, 2), p, DirectSol(1)), initialC);
q.compute();
%% VISUALIZE
v = Visualizer.create(fes);
prey = zeros(q.nT,1); pred = zeros(q.nT,1);
for k = 1:q.nT
  sol1 = q.history{k}(1:end/2);
  sol2 = q.history{k}(end/2+1:end);
  if ~mod(k,10)
    figure(1), clf, v.patch(sol1);
    axis([0 1 0 1 0 10]); caxis([0 2]), view(2)
    figure(2), clf, v.patch(sol2); colormap hot
    axis([0 1 0 1 0 10]); caxis([0 0.8]), view(2)
    fprintf('timestep: %d / %d\n', k, q.nT);
    drawnow;
  end
  prey(k) = mean(sol1); pred(k) = mean(sol2);
%  figure(3)
%  plot(1:k, prey(1:k), 1:k, pred(1:k)); axis([0 q.nT 0 1.5]);
%  title(sprintf('time %6.2f / %6.2f',(k-1)*T/q.nT, T));
end