classdef QpH1 < HierarchicElement
  methods % constructor
    function obj = QpH1(dim, order)
      obj = obj@HierarchicElement(dim, 2.^(1:dim), (order+1).^(1:dim), order);
      obj.doFTuple = (obj.order-1).^(0:dim);
      obj.conformity = 'H1';
    end
  end
  methods % evaluation
    function B = evalD0Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1));
      p = obj.order;
      points = 2*points-1; % transform to [-1,1]
      % shape functions
      N = cell(nP, 1);
      for d = 1:nP
        N{d} = cell(p+1,1);
        for i = 1:p+1
          N{d}{i} = obj.getShapeFunctions(points(:,d), i-1, 0);
        end
      end
      switch nP
        case 1
          for i = 1:p+1
            B(i,:) =  N{1}{i};
          end
        case 2
          offset = 0;
          % vertices
          for i = 1:2
            for j = 1:2
              B(offset+1,:) = N{1}{j}.*N{2}{i};
              offset = offset + 1;
            end
          end
          % faces
          for i = 0:1
            for j = 1:2
              for k = 1:p-1
                B(offset+1,:) =  N{mod(i,2)+1}{k+2}.*N{mod(i+1,2)+1}{j};
                offset = offset + 1;
              end
            end
          end
          % inner
          for j = 1:p-1
            for i = 1:p-1
              B(offset+1,:) =  N{1}{i+2}.*N{2}{j+2};
              offset = offset + 1;
            end
          end
        case 3
          offset = 0;
          % vertices
          for k = 1:2
            for j = 1:2
              for i = 1:2
                B(offset+1,:) = N{1}{i}.*N{2}{j}.*N{3}{k};
                offset = offset + 1;
              end
            end
          end
          % edges
          for i = 0:2
            for j = 1:2
              for k = 1:2
                for l = 1:p-1
                  B(offset+1,:) =  N{mod(i,3)+1}{l+2}.*N{mod(i+1,3)+1}{k}.*N{mod(i+2,3)+1}{j};
                  offset = offset + 1;
                end
              end
            end
          end
          % faces
          for i = 0:2
            for j = 1:2
              for k = 1:p-1
                for l = 1:p-1
                  B(offset+1,:) =  N{mod(i,3)+1}{l+2}.*N{mod(i+1,3)+1}{k+2}.*N{mod(i+2,3)+1}{j};
                  offset = offset + 1;
                end
              end
            end
          end
          % inner
          for k = 1:p-1
            for j = 1:p-1
              for i = 1:p-1
                B(offset+1,:) =  N{1}{i+2}.*N{2}{j+2}.*N{3}{k+2};
                offset = offset + 1;
              end
            end
          end
      end
    end
    function B = evalD1Basis(obj,points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP);
      p = obj.order;
      points = 2*points-1;
      % shape functions
      N = cell(nP, 1);
      dN = cell(nP, 1);
      for d = 1:nP
        N{1}{d} = cell(p+1,1);
        N{2}{d} = cell(p+1,1);
        for i = 1:p+1
          N{1}{d}{i} = obj.getShapeFunctions(points(:,d), i-1, 0);
          N{2}{d}{i} = 2*obj.getShapeFunctions(points(:,d), i-1, 1);
        end
      end
      switch nP
        case 1
          for i = 1:p+1
              B(i,:) =  N{2}{1}{i};
          end
        case 2
          offset = 0;
          % vertices
          for i = 1:2
            for j = 1:2
              B(offset+1,:,:) = [N{2}{1}{j}.*N{1}{2}{i}, N{1}{1}{j}.*N{2}{2}{i}];
              offset = offset + 1;
            end
          end
          % faces
          for i = 0:1
            for j = 1:2
              for k = 1:p-1
                B(offset+1,:,:) =  [N{(i==0) + 1}{mod(i,2)+1}{k+2}.*N{(i~=0) + 1}{mod(i+1,2)+1}{j}, ...
                                    N{(i==1) + 1}{mod(i,2)+1}{k+2}.*N{(i~=1) + 1}{mod(i+1,2)+1}{j}];
                offset = offset + 1;
              end
            end
          end
          % inner
          for j = 1:p-1
            for i = 1:p-1
              B(offset+1,:,:) =  [N{2}{1}{i+2}.*N{1}{2}{j+2}, N{1}{1}{i+2}.*N{2}{2}{j+2}];
              offset = offset + 1;
            end
          end
        case 3
          offset = 0;
          % vertices
          for k = 1:2
            for j = 1:2
              for i = 1:2
                B(offset+1,:,:) = [N{2}{1}{i}.*N{1}{2}{j}.*N{1}{3}{k}, ...
                                   N{1}{1}{i}.*N{2}{2}{j}.*N{1}{3}{k}, ...
                                   N{1}{1}{i}.*N{1}{2}{j}.*N{2}{3}{k}];
                offset = offset + 1;
              end
            end
          end
          % edges
          for i = 0:2
            for j = 1:2
              for k = 1:2
                for l = 1:p-1
                  B(offset+1,:,:) = [N{(i==0) + 1}{mod(i,3)+1}{l+2}.*N{(i==2) + 1}{mod(i+1,3)+1}{k}.*N{(i==1) + 1}{mod(i+2,3)+1}{j}, ...
                                     N{(i==1) + 1}{mod(i,3)+1}{l+2}.*N{(i==0) + 1}{mod(i+1,3)+1}{k}.*N{(i==2) + 1}{mod(i+2,3)+1}{j}, ...
                                     N{(i==2) + 1}{mod(i,3)+1}{l+2}.*N{(i==1) + 1}{mod(i+1,3)+1}{k}.*N{(i==0) + 1}{mod(i+2,3)+1}{j}];
                  offset = offset + 1;
                end
              end
            end
          end
          % faces
          for i = 0:2
            for j = 1:2
              for k = 1:p-1
                for l = 1:p-1
                  B(offset+1,:,:) = [N{(i==0) + 1}{mod(i,3)+1}{l+2}.*N{(i==2) + 1}{mod(i+1,3)+1}{k+2}.*N{(i==1) + 1}{mod(i+2,3)+1}{j}, ...
                                     N{(i==1) + 1}{mod(i,3)+1}{l+2}.*N{(i==0) + 1}{mod(i+1,3)+1}{k+2}.*N{(i==2) + 1}{mod(i+2,3)+1}{j}, ...
                                     N{(i==2) + 1}{mod(i,3)+1}{l+2}.*N{(i==1) + 1}{mod(i+1,3)+1}{k+2}.*N{(i==0) + 1}{mod(i+2,3)+1}{j}];
                  offset = offset + 1;
                end
              end
            end
          end
          % inner
          for k = 1:p-1
            for j = 1:p-1
              for i = 1:p-1
                B(offset+1,:,:) = [N{2}{1}{i+2}.*N{1}{2}{j+2}.*N{1}{3}{k+2}, ...
                                   N{1}{1}{i+2}.*N{2}{2}{j+2}.*N{1}{3}{k+2}, ...
                                   N{1}{1}{i+2}.*N{1}{2}{j+2}.*N{2}{3}{k+2}];
                offset = offset + 1;
              end
            end
          end
      end
      B = permute(B, [1 2 4 3]);
    end
    function B = evalD2Basis(obj,points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP, nP);
      p = obj.order;
      points = 2*points-1;
      switch nP
        case 1
          for i = 1:p+1
            B(i,:) =  4*obj.getShapeFunctions(2*points-1, i-1, 2);
          end
        case 2
          Btmp = zeros(obj.nB(nP), size(points,1), 3);
          % shape functions
          for d = 1:obj.dimension
            N{d} = cell(p+1,1);
            dN{d} = cell(p+1,1);
            d2N{d} = cell(p+1,1);
            for i = 1:p+1
              N{d}{i} = obj.getShapeFunctions(points(:,d), i-1, 0);
              dN{d}{i} = 2*obj.getShapeFunctions(points(:,d), i-1, 1);
              d2N{d}{i} = 4*obj.getShapeFunctions(points(:,d), i-1, 2);
            end
          end
          % vertices
          Btmp(1,:,3) = 1;
          Btmp(2,:,3) = -1;
          Btmp(3,:,3) = -1;
          Btmp(4,:,3) = 1;
          offset = 4;
          % faces
          for i = 0:1
            for j = 1:2
              for k = 1:p-1
                Btmp(offset+1,:,:) =  [d2N{mod(i,2)+1}{k+2}.*N{mod(i+1,2)+1}{j}, ...
                                       N{mod(i,2)+1}{k+2}.*d2N{mod(i+1,2)+1}{j}, ...
                                       dN{mod(i,2)+1}{k+2}.*dN{mod(i+1,2)+1}{j}];
                if i == 1
                  Btmp(offset+1,:,:) = Btmp(offset+1,:,[2 1 3]);
                end
                offset = offset + 1;
              end
            end
          end
          % inner
          for j = 1:p-1
            for i = 1:p-1
              Btmp(offset+1,:,:) =  [d2N{1}{i+2}.*N{2}{j+2}, N{1}{i+2}.*d2N{2}{j+2}, dN{1}{i+2}.*dN{2}{j+2}];
              offset = offset + 1;
            end
          end
          B(:,:,1,1) = Btmp(:,:,1);
          B(:,:,2,2) = Btmp(:,:,2);
          B(:,:,1,2) = Btmp(:,:,3);
          B(:,:,2,1) = Btmp(:,:,3);
        case 3
          error('! evalD2Basis: TODO !');
      end
      B = permute(B, [1 2 5 3 4]);
    end
  end
end