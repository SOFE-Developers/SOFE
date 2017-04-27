classdef GaussHex < QuadRule
  properties
  end
  methods % constructor
    function obj = GaussHex(n)
        obj = obj@QuadRule(n);
    end
  end
  methods % get methods
    function initData(obj)
      [p,w] = obj.getGaussPoints(); 
      p = 0.5 + 0.5*p;
      w = 0.5*w;
      wArray = bsxfun(@times,w*w',permute(w,[2 3 1]));
      obj.weights = wArray(:);
      obj.points = [kron(ones(length(w),1), kron(ones(length(w),1),p)) ...
                    kron(ones(length(w),1), kron(p, ones(length(w),1))) ...
                    kron(p, kron(ones(length(w),1), ones(length(w),1)))];
    end
  end
end