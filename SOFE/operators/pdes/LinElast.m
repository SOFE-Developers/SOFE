classdef LinElast < PDE
  methods % constructor
    function obj = LinElast(data, fes)
      mu = data.E/2/(1+data.nu);
      lambda = data.E*data.nu/(1+data.nu)/(1-2*data.nu);
      sLap = Op_data_SYMGRAD_SYMGRAD(2*mu, fes);
      divDiv = Op_data_div_div(lambda, fes);
      try
        f = Fc_Data_Id(data.f, fes, 0);
      catch
        f = Fc_Data_Id(@(x)0*x(:,1), fes, 0);
      end
      %
      lhs = {sLap, divDiv};
      rhs = {f};
      try % Neumann BC
        rhs{2} = Fc_Data_Id(data.g, fes, 1);
      end
      obj = obj@PDE({lhs}, {rhs});
    end
  end
end
