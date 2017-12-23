classdef FcId < Functional % ( F, V )
  methods % constructor
    function obj = FcId(coeff, fes, codim, varargin) % [loc]
      obj = obj@Functional(coeff, fes, codim, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      [points, weights] = obj.fes.getQuadData(obj.codim); % nPx1
      if isnumeric(obj.data)
        R = obj.fes.evalDoFVector(obj.data, [], obj.codim, 0, {k}); % nExnPx(nD*nD)
      else
        S = obj.observers{1}.evalState(k);
        R = obj.fes.evalFunction(obj.data, [], obj.codim, S, {k}); % nExnP
      end
      if ~isempty(points)
        basis = obj.fes.evalGlobalBasis([], obj.codim, 0, {k}); % nExnBxnPxnC
        [~,~,jac] = obj.fes.evalTrafoInfo([], obj.codim, {k}); % nExnP
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
        R = sum(bsxfun(@times, permute(R, [1 4 2 3]), basis), 4); % nExnBxnP
        R = sum(bsxfun(@times, R, permute(dX, [1 3 2])), 3); % nExnB
      end
    end
  end
end
