classdef CR < Element
  methods % constructor
    function obj = CR()
      obj = obj@Element(2, [2,3], [1,3], 1);
      obj.doFTuple = [0, 1, 0];
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
          B(1,:) = 1;
        case 2
          B(1,:) = 1-2*points(:,2);
          B(2,:) = 2*points(:,1)+2*points(:,2)-1;
          B(3,:) = 1-2*points(:,1);
      end
    end
    function B = evalD1Basis(obj,points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP);
      switch nP
        case 1
        case 2
          B(1,:,2) = -2;
          B(2,:,:) =  2;
          B(3,:,1) = -2;
      end
      B = permute(B, [1 2 4 3]);
    end
  end
end