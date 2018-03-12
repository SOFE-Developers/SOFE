classdef MassStokes < PDE
  methods % constructor
    function obj = MassStokes(C, FES, fes)
      opList = cell(3,1); lhs.sys = cell(3); rhs.sys = cell(3, 1);
      opList{1} = OpIdId(C, 0, FES, FES);
      opList{2} = OpIdId(0.0, 0, fes, fes);
      opList{3} = OpScalar(0.0);
      lhs.sys{1,1} = {1};
      lhs.sys{2,2} = {2};
      lhs.sys{3,3} = {3};
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end