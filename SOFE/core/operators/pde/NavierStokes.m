classdef NavierStokes < PDE
  methods % constructor
    function obj = NavierStokes(data, fesV, fesP)
%       lap = OpGradGrad(data.nu, fesV, fesV);
      lap = OpSGradSGrad(data.nu, fesV);
      conv = OpGradId(@(x,t,U)U{1}, fesV, fesV);
      grad = OpIdDiv(@(x)1+0*x(:,1), fesP, fesV);
      div = OpDivId(@(x)1+0*x(:,1), fesV, fesP);
      try
        f = {FcId(data.f, fesV, 0)};
      catch
        f = {};
      end
      %
      lhs = {{lap, conv},  {grad}; ...
             {div},     []};
      rhs = {f; []};
      obj = obj@PDE(lhs, rhs);
    end
  end
end
