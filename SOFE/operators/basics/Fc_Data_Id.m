classdef Fc_Data_Id < Functional % ( F, V )
  methods % constructor
    function obj = Fc_Data_Id(coeff, feSpace, codim, varargin) % [loc]
      obj = obj@Functional(coeff, feSpace, codim, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      [points, weights] = obj.feSpace.getQuadData(obj.codim); % nPx1
      if isreal(obj.data)
        R = obj.data;
      else
        R = obj.feSpace.evalFunction(obj.data, [], obj.codim, obj.state, k); % nExnPxnC
      end
      if ~isempty(points)
        basis = obj.feSpace.evalGlobalBasis([], obj.codim, 0, k); % nExnBxnPxnC
        [~,~,jac] = obj.feSpace.evalTrafoInfo([], obj.codim, k); % nExnP
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
        R = sum(bsxfun(@times, permute(R, [1 4 2 3]), basis), 4); % nExnBxnP
        R = sum(bsxfun(@times, R, permute(dX, [1 3 2])), 3); % nExnB
      end
    end
  end
end
