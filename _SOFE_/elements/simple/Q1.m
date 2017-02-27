classdef Q1 < Element
  methods % constructor
    function obj = Q1(dim)
      obj = obj@Element(dim, 2.^(1:dim), 2.^(1:dim), 1);
      obj.doFTuple = zeros(1,dim+1);
      obj.doFTuple(1) = 1;;
      obj.conformity = 'H1';
      obj.isLagrange = true;
    end
  end
  methods % evaluation
    function B = evalD0Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1));
      switch nP
        case 1
          B(1,:) = 1-points(:,1);
          B(2,:) = points(:,1);
        case 2
          B(1,:) = (1-points(:,1)).*(1-points(:,2));
          B(2,:) = points(:,1).*(1-points(:,2));
          B(3,:) = (1-points(:,1)).*points(:,2);
          B(4,:) = points(:,1).*points(:,2);
        case 3
          B(1,:) = (1-points(:,1)).*(1-points(:,2)).*(1-points(:,3));
          B(2,:) = points(:,1).*(1-points(:,2)).*(1-points(:,3));
          B(3,:) = (1-points(:,1)).*points(:,2).*(1-points(:,3));
          B(4,:) = points(:,1).*points(:,2).*(1-points(:,3));
          B(5,:) = (1-points(:,1)).*(1-points(:,2)).*points(:,3);
          B(6,:) = points(:,1).*(1-points(:,2)).*points(:,3);
          B(7,:) = (1-points(:,1)).*points(:,2).*points(:,3);
          B(8,:) = points(:,1).*points(:,2).*points(:,3);
      end
    end
    function B = evalD1Basis(obj, points)
      [nP, nD] = size(points);
      B = zeros(obj.nB(nD), nP, nD);
      switch nD
        case 1
          B(1,:,1) = -1;
          B(2,:,1) = 1;
        case 2
          B(1,:,1) = -(1-points(:,2));
          B(1,:,2) = -(1-points(:,1));
          B(2,:,1) = 1-points(:,2);
          B(2,:,2) = -points(:,1);
          B(3,:,1) = -points(:,2);
          B(3,:,2) = 1-points(:,1);
          B(4,:,1) = points(:,2);
          B(4,:,2) = points(:,1);
        case 3
          B(1,:,1) = -(1-points(:,2)).*(1-points(:,3));
          B(1,:,2) = -(1-points(:,1)).*(1-points(:,3));
          B(1,:,3) = -(1-points(:,1)).*(1-points(:,2));
          %
          B(2,:,1) = (1-points(:,2)).*(1-points(:,3));
          B(2,:,2) = -points(:,1).*(1-points(:,3));
          B(2,:,3) = -points(:,1).*(1-points(:,2));
          %
          B(3,:,1) = -points(:,2).*(1-points(:,3));
          B(3,:,2) = (1-points(:,1)).*(1-points(:,3));
          B(3,:,3) = -(1-points(:,1)).*points(:,2);
          %
          B(4,:,1) = points(:,2).*(1-points(:,3));
          B(4,:,2) = points(:,1).*(1-points(:,3));
          B(4,:,3) = -points(:,1).*points(:,2);
          %
          B(5,:,1) = -(1-points(:,2)).*points(:,3);
          B(5,:,2) = -(1-points(:,1)).*points(:,3);
          B(5,:,3) = (1-points(:,1)).*(1-points(:,2));
          %
          B(6,:,1) = (1-points(:,2)).*points(:,3);
          B(6,:,2) = -points(:,1).*points(:,3);
          B(6,:,3) = points(:,1).*(1-points(:,2));
          %
          B(7,:,1) = -points(:,2).*points(:,3);
          B(7,:,2) = (1-points(:,1)).*points(:,3);
          B(7,:,3) = (1-points(:,1)).*points(:,2);
          %
          B(8,:,1) = points(:,2).*points(:,3);
          B(8,:,2) = points(:,1).*points(:,3);
          B(8,:,3) = points(:,1).*points(:,2);
      end
      B = permute(B, [1 2 4 3]);
    end
  end
end