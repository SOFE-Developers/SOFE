classdef LinElast < PDE
  methods % constructor
    function obj = LinElast(data, fes)
      mu = data.E/2/(1+data.nu);
      lambda = data.E*data.nu/(1+data.nu)/(1-2*data.nu);
      sLap = OpSGradSGrad(2*mu, fes);
      divDiv = OpDivDiv(lambda, fes);
      try
        f = FcId(data.f, fes, 0);
      catch
        f = FcId(@(x)0*x(:,1), fes, 0);
      end
      %
      lhs = {sLap, divDiv};
      rhs = {f};
      try % Neumann BC
        rhs{2} = FcId(data.g, fes, 1);
      end
      obj = obj@PDE({lhs}, {rhs});
    end
  end
end
