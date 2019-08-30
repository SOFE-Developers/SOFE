classdef OpMatIdId < Operator % ( c*U, V )
  methods % constructor
    function obj = OpMatIdId(coeff, codim, fesTrial, varargin) % [fesTest loc]
      obj = obj@Operator(coeff, fesTrial, varargin{:});
      obj.codim = codim;
      if codim == 1 && isempty(obj.loc)
        obj.loc = @(x)~obj.fesTest.fixB(x);
      end
      obj.hasCoeff = false;
    end
  end
  methods
    function R = assembleOp(obj, k)
      basisJ = obj.fesTrial.evalGlobalBasis([], obj.codim, 0, {k}); % nExnBxnPxnW
      basisI = obj.fesTest.evalGlobalBasis([], obj.codim, 0, {k}); % nExnBxnPxnW
      %
      if isnumeric(obj.data)
        if size(obj.data,1)==1
          C = permute(obj.data, [3 1 2]); % 1x1xnC
        else
          C = obj.fesTrial.evalDoFVector(obj.data, [], 0, 0, {k}); % nExnPxnC
        end
        nW = size(obj.data);
      else
        S = obj.observers{1}.evalState(k);
        C = obj.fesTrial.evalFunction(obj.data, [], 0, S, {k}); % nExnPxnC
        nW = size(obj.data(zeros(1,obj.fesTrial.mesh.dimW)));
      end
      C = permute(C, [1 4 2 3]); % nEx1xnPxnC
      sz = size(C); sz = [sz(1:3) nW(:)'];
      C = reshape(C, sz); % nEx1xnPxnWxnW
      basisJ = sum(bsxfun(@times, permute(basisJ,[1 2 3 5 4]), C), 5); % nExnBxnPxnW
      %     
      R = obj.integrate(basisI, basisJ, k);
    end
  end
end