classdef FcCurl < Functional % ( F, curl(V) )
  methods % constructor
    function obj = FcCurl(coeff, fes, varargin) % [loc]
      obj = obj@Functional(coeff, fes, 0, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      [points, weights] = obj.fes.getQuadData(obj.codim); % nPx1
      if isnumeric(obj.data)
        R = obj.fes.evalDoFVector(obj.data, [], obj.codim, 0, {k}); % nExnPx(nD*nD)
      else
        try, S = obj.observers{1}.evalState(k); catch, S = obj.state; end
        R = obj.fes.evalFunction(obj.data, [], obj.codim, S, {k}); % nExnP
      end
      if ~isempty(points)
        basis = obj.fes.evalGlobalBasis([], obj.codim, 1, {k}); % nExnBxnPxnCxnD
        if obj.fes.element.dimension == 2
          curlBasis = basis(:,:,:,2,1) - basis(:,:,:,1,2);
        else
          curlBasis(:,:,:,3) = basis(:,:,:,2,1) - basis(:,:,:,1,2); % nExnBxnP
          curlBasis(:,:,:,2) = basis(:,:,:,1,3) - basis(:,:,:,3,1); % nExnBxnP
          curlBasis(:,:,:,1) = basis(:,:,:,3,2) - basis(:,:,:,2,3); % nExnBxnP
        end
        %
        [~,~,jac] = obj.fes.evalTrafoInfo([], obj.codim, {k}); % nExnP
        dX = bsxfun(@times, abs(jac), weights'); % nExnP
        R = sum(bsxfun(@times, permute(R, [1 4 2 3]), curlBasis), 4); % nExnBxnP
        R = sum(bsxfun(@times, R, permute(dX, [1 3 2])), 3); % nExnB
      end
    end
  end
end
