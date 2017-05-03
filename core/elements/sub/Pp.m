classdef Pp < Element
  methods % constructor
    function obj = Pp(dim, order)
      nB = [order+1, ...
           (order+1)*(order+2)/2, ...
           (order+1)*(order+2)*(order+3)/6];
      obj = obj@Element(dim, 2:(dim+1), nB(1:dim), order);
      obj.doFTuple(dim+1) = obj.nB(dim);
      obj.conformity = 'L2';
    end
  end
  methods % evaluation
    function B = evalD0Basis(obj, points)
      nD = size(points, 2);
      B = zeros(obj.nB(nD), size(points,1));
      p = obj.order;
      points = 2*points-1; % transform to [-1,1]^2
      N = cell(nD, 1);
      for d = 1:nD
        N{d} = cell(p+1,1);
        for i = 1:p+1
          N{d}{i} = obj.getLegendreFunctions(points(:,d), i-1, 0);
        end
      end
      switch nD
        case 1
          for i = 0:p
            B(i+1,:) =  N{1}{i+1};
          end
        case 2
          offset = 0;
          for deg = 0:p
            for i = 0:deg
              j = deg-i;
              B(offset+1,:) = N{1}{i+1}.*N{2}{j+1};
              offset = offset + 1;
            end
          end
        case 3
          offset = 0;
          for deg = 0:p
            for i = 0:deg
              for j = 0:deg-i
                k = deg-i-j;
                B(offset+1,:) = N{1}{i+1}.*N{2}{j+1}.*N{3}{k+1};
                offset = offset + 1;
              end
            end
          end
      end
    end
    function B = evalD1Basis(obj, points)
      nD = size(points, 2);
      B = zeros(obj.nB(nD), size(points,1), nD);
      p = obj.order;
      points = 2*points-1; % transform to [-1,1]^2
      N = cell(nD, 1);
      for d = 1:nD
        N{d} = cell(p+1,1);
        dN{d} = cell(p+1,1);
        for i = 1:p+1
          N{d}{i} = obj.getLegendreFunctions(points(:,d), i-1, 0);
          dN{d}{i} = 2*obj.getLegendreFunctions(points(:,d), i-1, 1);
        end
      end
      switch nD
        case 1
          for i = 0:p
            B(i+1,:) =  dN{1}{i+1};
          end
        case 2
          offset = 0;
          for deg = 0:p
            for i = 0:deg
              B(offset+1,:) = [dN{1}{i+1}.*N{2}{deg-i+1}; N{1}{i+1}.*dN{2}{deg-i+1}];
              offset = offset + 1;
            end
          end
        case 3
          offset = 0;
          for deg = 0:p
            for i = 0:deg
              for j = 0:deg-i
                k = deg-i-j;
                B(offset+1,:) = [dN{1}{i+1}.*N{2}{j+1}.*N{3}{k+1}; ...
                                 N{1}{i+1}.*dN{2}{j+1}.*N{3}{k+1}; ...
                                 N{1}{i+1}.*N{2}{j+1}.*dN{3}{k+1}];
                offset = offset + 1;
              end
            end
          end
      end
      B = permute(B, [1 2 4 3]);
    end
    function B = evalD2Basis(obj, points)
      nD = size(points, 2);
      B = zeros(obj.nB(nD), size(points,1), nD, nD);
      p = obj.order;
      points = 2*points-1; % transform to [-1,1]^2
      N = cell(nD, 1);
      for d = 1:nD
        N{d} = cell(p+1,1);
        dN{d} = cell(p+1,1);
        d2N{d} = cell(p+1,1);
        for i = 1:p+1
          N{d}{i} = obj.getLegendreFunctions(points(:,d), i-1, 0);
          dN{d}{i} = 2*obj.getLegendreFunctions(points(:,d), i-1, 1);
          d2N{d}{i} = 4*obj.getLegendreFunctions(points(:,d), i-1, 2);
        end
      end
      switch nD
        case 1
          for i = 0:p
            B(i+1,:) =  d2N{1}{i+1};
          end
        case 2
          offset = 0;
          for deg = 0:p
            for i = 0:deg
              B(offset+1,:,1,1) = d2N{1}{i+1}.*N{2}{deg-i+1};
              B(offset+1,:,2,2) = N{1}{i+1}.*d2N{2}{deg-i+1};
              B(offset+1,:,1,2) = dN{1}{i+1}.*dN{2}{deg-i+1};
              B(offset+1,:,2,1) = B(offset+1,1,2);
              offset = offset + 1;
            end
          end
        case 3
          offset = 0;
          for deg = 0:p
            for i = 0:deg
              for j = 0:deg-i
                k = deg-i-j;
                B(offset+1,:,1,1) = d2N{1}{i+1}.*N{2}{j+1}.*N{3}{k+1};
                B(offset+1,:,2,2) = N{1}{i+1}.*d2N{2}{j+1}.*N{3}{k+1};
                B(offset+1,:,3,3) = N{1}{i+1}.*N{2}{j+1}.*d2N{3}{k+1};
                B(offset+1,:,1,2) = dN{1}{i+1}.*dN{2}{j+1}.*N{3}{k+1};
                B(offset+1,:,1,3) = dN{1}{i+1}.*N{2}{j+1}.*dN{3}{k+1};
                B(offset+1,:,2,3) = N{1}{i+1}.*dN{2}{j+1}.*dN{3}{k+1};
                B(offset+1,:,2,1) = B(offset+1,1,2);
                B(offset+1,:,3,1) = B(offset+1,1,3);
                B(offset+1,:,3,2) = B(offset+1,2,3);
                offset = offset + 1;
              end
            end
          end
      end
      B = permute(B, [1 2 5 3 4]); % nExnPx1xnDxnD
    end
  end
end