classdef Qp < Element
  methods % constructor
    function obj = Qp(dim, order)
      if numel(order)==1, order = repmat(order,dim,1); end
      obj = obj@Element(dim, 2.^(1:dim), cumprod(order+1), order);
      obj.doFTuple(dim+1) = prod(obj.order+1);
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
        N{d} = cell(max(p)+1,1);
        for i = 1:max(p)+1
          N{d}{i} = obj.getLegendreFunctions(points(:,d), i-1, 0);
        end
      end
      switch nD
        case 1
          for i = 0:p
            B(i+1,:) =  N{1}{i+1};
          end
        case 2
          for j = 1:p(2)+1
            for i = 1:p(1)+1
              B(i+(p(1)+1)*(j-1),:) =  N{1}{i}.*N{2}{j};
            end
          end
        case 3
          for k = 1:p(3)+1
            for j = 1:p(2)+1
              for i = 1:p(1)+1
                B(i+(p(1)+1)*((j-1)+(p(2)+1)*(k-1)),:) =  N{1}{i}.*N{2}{j}.*N{3}{k};
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
        N{d} = cell(max(p)+1,1);
        dN{d} = cell(max(p)+1,1);
        for i = 1:max(p)+1
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
          for j = 1:p(2)+1
            for i = 1:p(1)+1
              B(i+(p(1)+1)*(j-1),:) =  [dN{1}{i}.*N{2}{j}; N{1}{i}.*dN{2}{j}];
            end
          end
        case 3
          for k = 1:p(3)+1
            for j = 1:p(2)+1
              for i = 1:p(1)+1
                B(i+(p(1)+1)*((j-1)+(p(2)+1)*(k-1)),:) =  [dN{1}{i}.*N{2}{j}.*N{3}{k}; ...
                                                           N{1}{i}.*dN{2}{j}.*N{3}{k}; ...
                                                           N{1}{i}.*N{2}{j}.*dN{3}{k}];
              end
            end
          end 
      end
      B = permute(B, [1 2 4 3]); % nExnPx1xnD
    end
  end
end