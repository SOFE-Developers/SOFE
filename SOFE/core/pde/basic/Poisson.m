classdef Poisson < PDE
  methods % constructor
    function obj = Poisson(data, fesTrial, varargin)
      if ~isempty(varargin), fesTest = varargin{:}; else, fesTest = fesTrial; end
      opList = {OpGradGrad(data.a, fesTrial, varargin{:}), ...
                FcId(data.f, fesTest, 0)};
      lhs.sys = {{1}};
      rhs.sys = {{2}};
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
  methods % residual
    function R = evalResidual(obj, U, points)
      assert(size(points,2)==obj.mesh.topology.dimP, '!Residual must be evaluated in elements!')
      f = obj.mesh.evalFunction(obj.list{2}.data, points, []); % nExnPxnC
      lapU = obj.fesTrial{1}.evalDoFVector(U, points, [], 2); % nExnPxnCx(nD*nD)
      iDiag = 1:(size(lapU,4)/2+1):size(lapU,4);
      R = f + sum(lapU(:,:,:,iDiag),4); % nExnPxnC
    end
  end
end