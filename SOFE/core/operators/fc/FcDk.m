classdef FcDk < Functional % ( F, dk(v) )
  properties
    k
  end
  methods % constructor
    function obj = FcDk(coeff, k, fes, varargin) % [loc]
      obj = obj@Functional(coeff, fes, 0, varargin{:});
      obj.k = k;
    end
  end
  methods
    function R = assembleOp(obj, k)
      [points, weights] = obj.fes.getQuadData(obj.codim); % nPx1     
      if isnumeric(obj.data)
        if size(obj.data,1)==1
          R = permute(obj.data, [3 1 2]); % 1x1xnC
        else
          R = obj.fes.evalDoFVector(obj.data, [], obj.codim, 0, {k}); % nExnPxnC
        end
      else
        try, S = obj.observers{1}.evalState(k); catch, S = obj.state; end
        R = obj.fes.evalFunction(obj.data, [], obj.codim, S, {k}); % nExnPxnC
      end
      if ~isempty(points)
        basis = obj.fes.evalGlobalBasis([], obj.codim, 1, {k}); % nExnBxnPx1xnD
        [~,~,jac] = obj.fes.evalTrafoInfo([], obj.codim, {k}); % nExnP
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
        R = sum(bsxfun(@times, permute(R, [1 4 2 3]), basis(:,:,:,:,obj.k)), 4); % nExnBxnP
        R = sum(bsxfun(@times, R, permute(dX, [1 3 2])), 3); % nExnB
      end
    end
  end
end
