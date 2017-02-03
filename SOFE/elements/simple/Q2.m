classdef Q2 < Element
  methods % constructor
      function obj = Q2()
          obj = obj@Element(2, [2,4], [3,9], 2);
          obj.doFTuple = [1, 1, 1];
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
          B(1,:) = (1-x).*(1-y);
          B(2,:) = x.*(1-y);
          B(3,:) = (1-x).*y;
          B(4,:) = x.*y;
          B(5,:) = 4*x.*(1-x).*(1-y);
          B(6,:) = 4*x.*y.*(1-x);
          B(7,:) = 4*(1-x).*y.*(1-y);
          B(8,:) = 4*x.*y.*(1-y);
          B(9,:) = 16*x.*(1-x).*y.*(1-y);
      end
    end
    function B = evalD1Basis(obj,points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP);
      switch nP
        case 1
          B(1,:,1) = -1;
          B(2,:,1) = 1;
          B(3,:,1) = 4*(1-2*points);
        case 2
          x = points(:,1); y = points(:,2);
          B(1,:,1) = -(1-y);
          B(1,:,2) = -(1-x);
          B(2,:,1) = 1-y;
          B(2,:,2) = -x;
          B(3,:,1) = -y;
          B(3,:,2) = 1-x;
          B(4,:,1) = y;
          B(4,:,2) = x;
          B(5,:,1) = 4*(1-2*x).*(1-y);
          B(5,:,2) = -4*x.*(1-x);
          B(6,:,1) = 4*y.*(1-2*x);
          B(6,:,2) = 4*x.*(1-x);
          B(7,:,1) = -4*y.*(1-y);
          B(7,:,2) = 4*(1-x).*(1-2*y);
          B(8,:,1) = 4*y.*(1-y);
          B(8,:,2) = 4*x.*(1-2*y);
          B(9,:,1) = 16*(1-2*x).*y.*(1-y);
          B(9,:,2) = 16*x.*(1-x).*(1-2*y);
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
          B(5,:,1,2) = -4*(1-2*y);
          B(5,:,2,2) = -8*(1-x);
          B(6,:,1,2) = 4*(1-2*y);
          B(6,:,2,2) = -8*x;
          B(7,:,1,1) = -8*(1-y);
          B(7,:,1,2) = -4*(1-2*x);
          B(8,:,1,1) = -8*y;
          B(8,:,1,2) = 4*(1-2*x);
          B(9,:,1,1) = -32*y.*(1-y);
          B(9,:,1,2) = 16*(1-2*x).*(1-2*y);
          B(9,:,2,2) = -32*x.*(1-x);
          %
          B(:,:,2,1) = B(:,:,1,2);
      end
      B = permute(B, [1 2 5 3 4]);
    end
  end
end