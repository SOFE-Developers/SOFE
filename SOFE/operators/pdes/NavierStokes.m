classdef NavierStokes < PDE
  methods % constructor
    function obj = NavierStokes(data, fesV, fesP)
      lap = Op_data_GRAD_GRAD(data.nu, fesV, fesV);
      conv = Op_Data_GRAD_Id(@(x,t,U)U{1}, fesV, fesV);
      grad = Op_data_Grad_Id(@(x)1+0*x(:,1), fesP, fesV);
      div = Op_data_div_id(@(x)1+0*x(:,1), fesV, fesP);
      try
        f = {Fc_Data_Id(data.f, fesV, 0)};
      catch
        f = {};
      end
      %
      lhs = {{lap, conv},  {grad}; ...
             {div},     []};
      rhs = {f; []};
      obj = obj@PDE(lhs, rhs, struct('nonLin', true));
    end
    function assemble(obj)
      assemble@PDE(obj);
      % fix pressure doF
      nDoFV = obj.fesTrial{1}.getNDoF();
      nDoFP = obj.fesTrial{2}.getNDoF();
      obj.fDoFsTrial(nDoFV+nDoFP) = false;
      obj.fDoFsTest(nDoFV+nDoFP) = false;
    end
  end
end
