classdef GaussQuad < QuadRule
  properties
  end
  methods % constructor
    function obj = GaussQuad(n,varargin)
        obj = obj@QuadRule(n);
    end
  end
  methods % get methods
    function initData(obj)
      [p,w] = obj.getGaussPoints(); 
      p = 0.5 + 0.5*p;
      w = 0.5*w;
      wArray = w*w';
      obj.weights = wArray(:);
      obj.points = [kron(ones(length(w),1),p) kron(p, ones(length(w),1))];
    end
  end
end