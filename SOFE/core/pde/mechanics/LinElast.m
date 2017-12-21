classdef LinElast < PDE
  methods % constructor
    function obj = LinElast(data, fes)
      mu = data.E/2/(1+data.nu);
      lambda = data.E*data.nu/(1+data.nu)/(1-2*data.nu);
      opList = {OpSGradSGrad(2*mu, fes), ...
                OpDivDiv(lambda, fes),[], []};
      try
        fData = data.f;
      catch
        fData = @(x)0*x(:,1);
      end
      opList{3} = FcId(fData, fes, 0);
      %
      lhs.sys = {{1, 2}};
      %
      rhs.sys = {{3}};
      try
        opList{4} = FcId(data.g, fes, 1);
        rhs.sys{1} = [rhs.sys{1} 4];
      catch
        opList = opList(1:3);
      end
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end
