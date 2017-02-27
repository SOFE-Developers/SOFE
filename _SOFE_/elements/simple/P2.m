classdef P2 < Element
  properties
  end
  methods % constructor
    function obj = P2()
      obj = obj@Element(2, [2,3], [3,6], 2);
      obj.doFTuple = [1 1 0];
      obj.conformity = 'H1';
      obj.isLagrange = false;
    end
  end
  methods % evaluation
    function B = evalD0Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1));
      switch nP
        case 1
          B(1,:) = 1-points;
          B(2,:) = points;
          B(3,:) = 4*points.*(1-points);
        case 2
          x = points(:,1); y = points(:,2);
          B(1,:) = 1-x-y;
          B(2,:) = x;
          B(3,:) = y;
          B(4,:) = 4*x.*(1-x-y);
          B(5,:) = 4*x.*y;
          B(6,:) = 4*y.*(1-x-y);
        case 3
      end
    end
    function B = evalD1Basis(obj,points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP);
      switch nP
        case 1
          B(1,:,1) = -1;
          B(2,:,1) = 1;
          B(3,:,1) = 4-8*points;
        case 2
          x = points(:,1); y = points(:,2);
          B(1,:,:) = -1;
          B(2,:,1) = 1;
          B(3,:,2) = 1;
          B(4,:,1) = 4*(1-2*x-y);
          B(4,:,2) = -4*x;
          B(5,:,1) = 4*y;
          B(5,:,2) = 4*x;
          B(6,:,1) = -4*y;
          B(6,:,2) = 4*(1-x-2*y);
        case 3
      end
      B = permute(B, [1 2 4 3]);
    end
    function B = evalD2Basis(obj,points)
      nP = size(points,2);
      B = zeros(obj.nB(nP), size(points,1), nP, nP);
      switch nP
        case 1
          B(3,:,1) = -8;
        case 2
          x = points(:,1); y = points(:,2);
          B(4,:,1,1) = -8;
          B(4,:,1,2) = -4;
          B(4,:,2,1) = -4;
          B(5,:,1,2) = 4;
          B(5,:,2,1) = 4;
          B(6,:,1,2) = -4;
          B(6,:,2,1) = -4;
          B(6,:,2,2) = -8;
        case 3
      end
      B = permute(B, [1 2 5 3 4]);
    end
  end
end