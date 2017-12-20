classdef CDR2 < PDE2
% GradGrad(a) + Grad(b) + Id(c) + Id_Gamma(h) = FId(f) + FId_Gamma(g)
  methods % constructor
    function obj = CDR2(data, fesTrial, varargin)
      if ~isempty(varargin), fesTest = varargin{:}; else, fesTest = fesTrial; end
      opList = cell(1,6);
      lhs.sys = {{}}; rhs.sys = {{}};
      %
      opList{1} = OpGradGrad(data.a, fesTrial, fesTest); nOp = 1;
      lhs.sys{1} = [lhs.sys{1} nOp];
      try %#ok<*TRYNC>
        opList{nOp + 1} = OpGradId(data.b, fesTrial, varargin{:}); nOp = nOp + 1;
        lhs.sys{1} = [lhs.sys{1} nOp];
      end
      try
        opList{nOp + 1} = OpIdId(data.c, 0, fesTrial, varargin{:}); nOp = nOp + 1;
        lhs.sys{1} = [lhs.sys{1} nOp];
      end
      try
        opList{nOp + 1} = OpIdId(data.h, 1, fesTrial, fesTest); nOp = nOp + 1;
        lhs.sys{1} = [lhs.sys{1} nOp];
      end
      %
      try fData = data.f; catch, fData = @(x)0*x(:,1);end
      opList{nOp + 1} = FcId(fData, fesTest, 0); nOp = nOp + 1;
      rhs.sys{1} = [rhs.sys{1} nOp];
      try
        opList{nOp + 1} = FcId(data.g, fesTest, 1); nOp = nOp + 1;
        rhs.sys{1} = [rhs.sys{1} nOp];
      end
      opList = opList(1:nOp);
      %
      obj = obj@PDE2(opList, lhs, rhs);
    end
  end
end
