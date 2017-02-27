classdef Op_DATA_GRAD_GRAD < Operator % ( CC*GRAD(U), GRAD(V) )
  methods % constructor
    function obj = Op_DATA_GRAD_GRAD(coeff, feSpaceTrial, varargin)
      obj = obj@Operator(coeff, feSpaceTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      gradBasisJ = obj.feSpaceTrial.evalGlobalBasis([], 0, 1, k); % nExnBxnPxnCxnD
      if obj.isGalerkin
        gradBasisI = gradBasisJ;
      else
        gradBasisI = obj.feSpaceTest.evalGlobalBasis([], 0, 1, k); % nExnBxnPxnCxnD
      end
      %
      [~, weights] = obj.feSpaceTrial.getQuadData(obj.codim);
      [~,~,jac] = obj.feSpaceTrial.evalTrafoInfo([], obj.codim, k);
      coeff = obj.feSpaceTrial.evalFunction(obj.data, [], obj.codim, obj.state, k);
      nD = obj.feSpaceTrial.element.dimension;
      coeff = reshape(coeff, size(coeff,1), size(coeff,2), nD, nD); % nExnPxnDxnD
      dX = bsxfun(@times, abs(jac), weights'); % nExnP
      nE = size(gradBasisI, 1); nBI = size(gradBasisI, 2);
      nBJ = size(gradBasisJ,2); nP = size(gradBasisI,3);
      %
      gradBasisJ = sum(bsxfun(@times, permute(coeff, [1 6 2 3 5 4]), ...
                                      permute(gradBasisJ, [1 2 3 6 4 5])), 6); % nExnBxnPxnDxnD
      %
      gradBasisI = reshape(gradBasisI, nE, nBI, nP, []);
      gradBasisJ = reshape(gradBasisJ, nE, nBJ, nP, []);
      R = sum(bsxfun(@times, permute(gradBasisI, [1 2 5 3 4]), ...
                             permute(gradBasisJ, [1 5 2 3 4])), 5); % nExnBIxnBJxnPx[]
      R = sum(bsxfun(@times, R, permute(dX, [1 3 4 2])), 4); % nExnBIxnBJ
    end
  end
end
