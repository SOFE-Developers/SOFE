classdef OpGradGrad < Operator % ( c*GRAD(U), GRAD(V) )
  methods % constructor
    function obj = OpGradGrad(coeff, fesTrial, varargin)
      obj = obj@Operator(coeff, fesTrial, varargin{:});
    end
  end
  methods
    function R = assembleOp(obj, k)
      gradBasisJ = obj.fesTrial.evalGlobalBasis([], 0, 1, {k}); % nExnBxnPxnCxnD
      gradBasisI = obj.fesTest.evalGlobalBasis([], 0, 1, {k});
      R = obj.integrate(true, gradBasisI, gradBasisJ, k);
    end
    function R = assembleRef_(obj) % deprecated
      [~, weights] = obj.fesTrial.getQuadData(obj.codim);
      [~,~,jac] = obj.fesTrial.evalTrafoInfo([], obj.codim, 1); % nExnP
      dX = bsxfun(@times, abs(jac), weights'); % nExnP
      basisI = obj.fesTrial.evalGlobalBasis([], 0, 1, 1); % nExnBxnPxnCxnD
      basisJ = obj.fesTest.evalGlobalBasis([], 0, 1, 1);
      nE = size(basisI, 1); nBI = size(basisI, 2); nBJ = size(basisJ,2); nP = size(basisI,3);
      basisI = reshape(basisI, nE, nBI, nP, []); basisJ = reshape(basisJ, nE, nBJ, nP, []);
      R = sum(bsxfun(@times, permute(basisI, [1 2 5 3 4]), ...
                             permute(basisJ, [1 5 2 3 4])), 5); % nExnBIxnBJxnPx(nC*nD)
      R = sum(bsxfun(@times, R, permute(dX, [1 3 4 2])), 4); % nExnBIxnBJ
      R = permute(R, [2 3 1]);
      R(abs(R)<1e-15) = 0;
    end
  end
end
