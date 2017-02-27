classdef P1 < Element
  methods % constructor
    function obj = P1(dim)
      obj = obj@Element(dim, 2:(dim+1), 2:(dim+1), 1);
      obj.doFTuple = zeros(1,dim+1);
      obj.doFTuple(1) = 1;
      obj.conformity = 'H1';
      obj.isLagrange = true;
    end
  end
  methods % local evaluation
    function B = evalD0Basis(obj, points)
      [nP, nD] = size(points);
      B = zeros(obj.nB(nD), nP);
      B(1,:) = 1-sum(points,2);
      for i = 1:nD
        B(i+1,:) = points(:,i);
      end
    end
    function B = evalD1Basis(obj, points)
      [nP, nD] = size(points);
      B = zeros(obj.nB(nD), nP, 1, nD);
      B(1,:,1,:) = -1;
      for i = 1:nD
        B(i+1,:,1,i) = 1;
      end
    end
  end
end