classdef PpH1 < HierarchicElement
  methods % constructor
    function obj = PpH1(dim, order)
      nB = [order+1, ...
           (order+1)*(order+2)/2, ...
           (order+1)*(order+2)*(order+3)/6 + 4*2*(order-1)*(order-2)/2];
      obj = obj@HierarchicElement(dim, 2:(dim+1), nB(1:dim), order);
      for i = 0:dim
        obj.doFTuple(i+1) = prod(order - (1:i))/factorial(i);
      end
      obj.conformity = 'H1';
    end
  end
  methods % local evaluation
    function B = evalD0Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1));
      p = obj.order;
      switch nP
        case 1
          for i = 1:p+1
            B(i,:) =  obj.getShapeFunctions(2*points-1, i-1, 0);
          end
        case 2
          L{1} = 1 - sum(points,2);
          L{2} = points(:,1);
          L{3} = points(:,2);
          % Kernel functions
          kernel = cell(3,1);
          for k = 1:3
            kernel{k} = cell(3,1);
          end
          kernel{2}{1} = obj.getKernel(L{2}-L{1});
          kernel{3}{2} = obj.getKernel(L{3}-L{2});
          kernel{1}{3} = obj.getKernel(L{3}-L{1});
          % vertices
          for i = 1:3
            B(i,:) = L{i};
          end         
          offset = 3;
          % faces
          for f = 1:3
            for i = 1:p-1
              B(offset+1,:) =  L{f}.*L{mod(f,3)+1}.*kernel{mod(f,3)+1}{f}(:,i);
              offset = offset + 1;
            end
          end
          % inner
          for j = 1:p-2
            for i = 1:p-1-j
              B(offset+1,:) =  L{1}.*L{2}.*L{3}.*kernel{2}{1}(:,i).*kernel{1}{3}(:,j);
              offset = offset + 1;
            end
          end
        case 3
          L1 = 1 - sum(points,2);
          L2 = points(:,1);
          L3 = points(:,2);
          L4 = points(:,3);
          % Kernel functions
          kernel2m1 = obj.getKernel(L2-L1);
          kernel3m2 = obj.getKernel(L3-L2);
          kernel3m1 = obj.getKernel(L3-L1);
          kernel4m1 = obj.getKernel(L4-L1);
          kernel4m2 = obj.getKernel(L4-L2);
          kernel4m3 = obj.getKernel(L4-L3);
          % vertices
          B(1,:) = L1;
          B(2,:) = L2;
          B(3,:) = L3;
          B(4,:) = L4;
          offset = 4;
          % edge
          for i = 1:p-1
              B(offset+1,:) =  L1.*L2.*kernel2m1(:,i);
              offset = offset + 1;
          end
          for i = 1:p-1
              B(offset+1,:) =  L2.*L3.*kernel3m2(:,i);
              offset = offset + 1;
          end
          for i = 1:p-1
              B(offset+1,:) =  L3.*L1.*kernel3m1(:,i);
              offset = offset + 1;
          end
          for i = 1:p-1
              B(offset+1,:) =  L1.*L4.*kernel4m1(:,i);
              offset = offset + 1;
          end
          for i = 1:p-1
              B(offset+1,:) =  L2.*L4.*kernel4m2(:,i);
              offset = offset + 1;
          end
          for i = 1:p-1
              B(offset+1,:) =  L3.*L4.*kernel4m3(:,i);
              offset = offset + 1;
          end
          % face
%           tmp = (p-1)*(p-2)/2;
          L123 = L1.*L2.*L3;
          L124 = L1.*L2.*L4;
          L234 = L2.*L3.*L4;
          L134 = L1.*L3.*L4;
%          for dg = 3:p
%              for i = 1:dg-2
          for j = 1:p-2
            for i = 1:p-1-j
                %j = dg-1-i;
                si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
                % origin @ vertex 1
                B(offset + 1,:) = L123.*kernel2m1(:,i).*kernel3m1(:,j);
                % origin @ vertex 2
                B(offset + 2,:) = sj*L123.*kernel3m2(:,i).*kernel2m1(:,j);
                % origin @ vertex 3
                B(offset + 3,:) = sij*L123.*kernel3m1(:,i).*kernel3m2(:,j);
                offset = offset + 3;
            end
          end
          for j = 1:p-2
            for i = 1:p-1-j
              si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
              % origin @ vertex 1
              B(offset + 1,:) = L124.*kernel2m1(:,i).*kernel4m1(:,j);
              % origin @ vertex 2
              B(offset + 2,:) = sj*L124.*kernel4m2(:,i).*kernel2m1(:,j);
              % origin @ vertex 3
              B(offset + 3,:) = sij*L124.*kernel4m1(:,i).*kernel4m2(:,j);
              offset = offset + 3;
            end
          end
          for j = 1:p-2
            for i = 1:p-1-j
              si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
              % origin @ vertex 1
              B(offset + 1,:) = L234.*kernel3m2(:,i).*kernel4m2(:,j);
              % origin @ vertex 2
              B(offset + 2,:) = sj*L234.*kernel4m3(:,i).*kernel3m2(:,j);
              % origin @ vertex 3
              B(offset + 3,:) = sij*L234.*kernel4m2(:,i).*kernel4m3(:,j);
              offset = offset + 3;
            end
          end
          for j = 1:p-2
            for i = 1:p-1-j
              si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
              % origin @ vertex 1
              B(offset + 1,:) = L134.*kernel3m1(:,i).*kernel4m1(:,j);
              % origin @ vertex 2
              B(offset + 2,:) = sj*L134.*kernel4m3(:,i).*kernel3m1(:,j);
              % origin @ vertex 3
              B(offset + 3,:) = sij*L134.*kernel4m1(:,i).*kernel4m3(:,j);
              offset = offset + 3;
            end
          end
          % inner
          L1234 = L123.*L4;
          for k = 1:p-3
            for j = 1:p-2-k
              for i = 1:p-1-j-k
                B(offset+1,:) =  L1234.*kernel2m1(:,i).*kernel3m1(:,j).*kernel4m1(:,k);
                offset = offset + 1;
              end
            end
          end
      end
    end
    function B = evalD1Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP);
      p = obj.order;
      switch nP
        case 1
          for i = 1:p+1
            B(i,:) =  2*obj.getShapeFunctions(2*points-1, i-1, 1);
          end
        case 2
          L{1} = 1 - sum(points,2);
          L{2} = points(:,1);
          L{3} = points(:,2);
          % Kernel functions
          kernel = cell(3,1); dKernel = cell(3,1);
          for k = 1:3
            kernel{k} = cell(3,1);
            dKernel{k} = cell(3,1);
          end
          [kernel{2}{1},dKernel{2}{1}] = obj.getKernel(L{2}-L{1});
          [kernel{3}{2},dKernel{3}{2}] = obj.getKernel(L{3}-L{2});
          [kernel{1}{3},dKernel{1}{3}] = obj.getKernel(L{3}-L{1});
          % vertices
          B(1,:,1) = -1;
          B(1,:,2) = -1;
          B(2,:,1) = 1;
          B(3,:,2) = 1;
          offset = 3;
          % faces          
          for i = 1:p-1
            B(offset+1,:,:) =  [kernel{2}{1}(:,i).*(L{1}-L{2})+2*L{1}.*L{2}.*dKernel{2}{1}(:,i), ...
                                   L{2}.*(-kernel{2}{1}(:,i)+L{1}.*dKernel{2}{1}(:,i))];
            offset = offset + 1;
          end
          for i = 1:p-1
            B(offset+1,:,:) =  [L{3}.*(kernel{3}{2}(:,i)-L{2}.*dKernel{3}{2}(:,i)), ...
                                    L{2}.*(kernel{3}{2}(:,i)+L{3}.*dKernel{3}{2}(:,i))];
            offset = offset + 1;
          end
          for i = 1:p-1
            B(offset+1,:,:) =  [L{3}.*(-kernel{1}{3}(:,i)+L{1}.*dKernel{1}{3}(:,i)), ...
                                   kernel{1}{3}(:,i).*(L{1}-L{3})+2*L{1}.*L{3}.*dKernel{1}{3}(:,i)];
            offset = offset + 1;
          end
          % inner
          for j = 1:p-2
            for i = 1:p-1-j
              B(offset+1,:,:) = ...
                [L{3}.*((L{1}-L{2}).*kernel{2}{1}(:,i).*kernel{1}{3}(:,j) + ...
                 L{1}.*L{2}.*(2*dKernel{2}{1}(:,i).*kernel{1}{3}(:,j) + kernel{2}{1}(:,i).*dKernel{1}{3}(:,j))), ...
                 L{2}.*((L{1}-L{3}).*kernel{2}{1}(:,i).*kernel{1}{3}(:,j) + ...
                 L{1}.*L{3}.*(dKernel{2}{1}(:,i).*kernel{1}{3}(:,j) + 2*kernel{2}{1}(:,i).*dKernel{1}{3}(:,j)))];
              offset = offset + 1;
            end
          end
        case 3
          L1 = 1 - sum(points,2);
          L2 = points(:,1);
          L3 = points(:,2);
          L4 = points(:,3);
          % Kernel functions
          [kernel2m1, dKernel2m1] = obj.getKernel(L2-L1);
          [kernel3m2, dKernel3m2] = obj.getKernel(L3-L2);
          [kernel3m1, dKernel3m1] = obj.getKernel(L3-L1);
          [kernel4m1, dKernel4m1] = obj.getKernel(L4-L1);
          [kernel4m2, dKernel4m2] = obj.getKernel(L4-L2);
          [kernel4m3, dKernel4m3] = obj.getKernel(L4-L3);
          % vertices
          B(1,:,1) = -1; B(1,:,2) = -1; B(1,:,3) = -1;
          B(2,:,1) = 1;
          B(3,:,2) = 1;
          B(4,:,3) = 1;
          offset =  4;
          % edge
          for i = 1:p-1
            B(offset+1,:,:) =  [(L1-L2).*kernel2m1(:,i) + 2*L1.*L2.*dKernel2m1(:,i), ...
                                   L2.*(-kernel2m1(:,i) + L1.*dKernel2m1(:,i)), ...
                                   L2.*(-kernel2m1(:,i) + L1.*dKernel2m1(:,i))];
            offset = offset + 1;
          end
          for i = 1:p-1
            B(offset+1,:,:) =  [L3.*(kernel3m2(:,i) - L2.*dKernel3m2(:,i)), ...
                                   L2.*(kernel3m2(:,i) + L3.*dKernel3m2(:,i)), ...
                                   0.*kernel3m2(:,i)];
            offset = offset + 1;
          end
          for i = 1:p-1
            B(offset+1,:,:) =  [-L3.*(kernel3m1(:,i) - L1.*dKernel3m1(:,i)), ...
                                   (L1-L3).*kernel3m1(:,i) + 2*L1.*L3.*dKernel3m1(:,i), ...
                                   -L3.*(kernel3m1(:,i) - L1.*dKernel3m1(:,i))];
            offset = offset + 1;
          end
          for i = 1:p-1
            B(offset+1,:,:) =  [L4.*(-kernel4m1(:,i) + L1.*dKernel4m1(:,i)), ...
                                   L4.*(-kernel4m1(:,i) + L1.*dKernel4m1(:,i)), ...
                                   (L1-L4).*kernel4m1(:,i) + 2*L1.*L4.*dKernel4m1(:,i)];
            offset = offset + 1;
          end
          for i = 1:p-1
            B(offset+1,:,:) =  [L4.*(kernel4m2(:,i) - L2.*dKernel4m2(:,i)), ...
                                   0.*kernel4m2(:,i), ...
                                   L2.*(kernel4m2(:,i) + L4.*dKernel4m2(:,i))];
            offset = offset + 1;
          end
          for i = 1:p-1
            B(offset+1,:,:) =  [0.*kernel4m3(:,i), ...
                                   L4.*(kernel4m3(:,i) - L3.*dKernel4m3(:,i)), ...
                                   L3.*(kernel4m3(:,i) + L4.*dKernel4m3(:,i))];
            offset = offset + 1;
          end
          % face
%           tmp = (p-1)*(p-2)/2;
          L123 = L1.*L2.*L3;
          L124 = L1.*L2.*L4;
          L234 = L2.*L3.*L4;
          L134 = L1.*L3.*L4;
          for j = 1:p-2
            for i = 1:p-1-j
              si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
              % ------------
              % origin @ vertex 1
              % ------------
              B(offset+1,:,1) = L3.*(L1-L2).*kernel2m1(:,i).*kernel3m1(:,j) + ...
                      L123.*(2*dKernel2m1(:,i).*kernel3m1(:,j) + kernel2m1(:,i).*dKernel3m1(:,j));
              B(offset+1,:,2) = L2.*(L1-L3).*kernel2m1(:,i).*kernel3m1(:,j) + ...
                      L123.*(dKernel2m1(:,i).*kernel3m1(:,j) + 2*kernel2m1(:,i).*dKernel3m1(:,j));
              B(offset+1,:,3) = -L2.*L3.*kernel2m1(:,i).*kernel3m1(:,j) + ...
                      L123.*(dKernel2m1(:,i).*kernel3m1(:,j) + kernel2m1(:,i).*dKernel3m1(:,j));
              % ------------
              % origin @ vertex 2
              % ------------          
              B(offset+2,:,1) = sj*(L3.*(L1-L2).*kernel3m2(:,i).*kernel2m1(:,j) + ...
                      L123.*(-dKernel3m2(:,i).*kernel2m1(:,j) + 2*kernel3m2(:,i).*dKernel2m1(:,j)));
              B(offset+2,:,2) = sj*(L2.*(L1-L3).*kernel3m2(:,i).*kernel2m1(:,j) + ...
                      L123.*(dKernel3m2(:,i).*kernel2m1(:,j) + kernel3m2(:,i).*dKernel2m1(:,j)));
              B(offset+2,:,3) = sj*(-L2.*L3.*kernel3m2(:,i).*kernel2m1(:,j) + ...
                      L123.*(0*dKernel3m2(:,i).*kernel2m1(:,j) + kernel3m2(:,i).*dKernel2m1(:,j)));
              % ------------
              % origin @ vertex 3
              % ------------      
              B(offset+3,:,1) = sij*(L3.*(L1-L2).*kernel3m1(:,i).*kernel3m2(:,j) + ...
                      L123.*(dKernel3m1(:,i).*kernel3m2(:,j) - kernel3m1(:,i).*dKernel3m2(:,j)));
              B(offset+3,:,2) = sij*(L2.*(L1-L3).*kernel3m1(:,i).*kernel3m2(:,j) + ...
                      L123.*(2*dKernel3m1(:,i).*kernel3m2(:,j) + kernel3m1(:,i).*dKernel3m2(:,j)));
              B(offset+3,:,3) = sij*(-L2.*L3.*kernel3m1(:,i).*kernel3m2(:,j) + ...
                      L123.*(dKernel3m1(:,i).*kernel3m2(:,j) + 0*kernel3m1(:,i).*dKernel3m2(:,j)));
              offset = offset + 3;
            end
          end
          % face 2
          for j = 1:p-2
            for i = 1:p-1-j
              si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
              % ------------
              % origin @ vertex 1
              % ------------
              B(offset+1,:,1) = L4.*(L1-L2).*kernel2m1(:,i).*kernel4m1(:,j) + ...
                      L124.*(2*dKernel2m1(:,i).*kernel4m1(:,j) + kernel2m1(:,i).*dKernel4m1(:,j));
              B(offset+1,:,2) = -L2.*L4.*kernel2m1(:,i).*kernel4m1(:,j) + ...
                      L124.*(dKernel2m1(:,i).*kernel4m1(:,j) + kernel2m1(:,i).*dKernel4m1(:,j));
              B(offset+1,:,3) = L2.*(L1-L4).*kernel2m1(:,i).*kernel4m1(:,j) + ...
                      L124.*(dKernel2m1(:,i).*kernel4m1(:,j) + 2*kernel2m1(:,i).*dKernel4m1(:,j));
              % ------------
              % origin @ vertex 2
              % ------------          
              B(offset+2,:,1) = sj*(L4.*(L1-L2).*kernel4m2(:,i).*kernel2m1(:,j) + ...
                      L124.*(-dKernel4m2(:,i).*kernel2m1(:,j) + 2*kernel4m2(:,i).*dKernel2m1(:,j)));
              B(offset+2,:,2) = sj*(-L2.*L4.*kernel4m2(:,i).*kernel2m1(:,j) + ...
                      L124.*(0*dKernel4m2(:,i).*kernel2m1(:,j) + kernel4m2(:,i).*dKernel2m1(:,j)));
              B(offset+2,:,3) = sj*(L2.*(L1-L4).*kernel4m2(:,i).*kernel2m1(:,j) + ...
                      L124.*(dKernel4m2(:,i).*kernel2m1(:,j) + kernel4m2(:,i).*dKernel2m1(:,j)));
              % ------------
              % origin @ vertex 3
              % ------------      
              B(offset+3,:,1) = sij*(L4.*(L1-L2).*kernel4m1(:,i).*kernel4m2(:,j) + ...
                      L124.*(dKernel4m1(:,i).*kernel4m2(:,j) - kernel4m1(:,i).*dKernel4m2(:,j)));
              B(offset+3,:,2) = sij*(-L2.*L4.*kernel4m1(:,i).*kernel4m2(:,j) + ...
                      L124.*(dKernel4m1(:,i).*kernel4m2(:,j) + 0*kernel4m1(:,i).*dKernel4m2(:,j)));
              B(offset+3,:,3) = sij*(L2.*(L1-L4).*kernel4m1(:,i).*kernel4m2(:,j) + ...
                      L124.*(2*dKernel4m1(:,i).*kernel4m2(:,j) + kernel4m1(:,i).*dKernel4m2(:,j)));
              offset = offset + 3;
            end
          end
          % face 3
          for j = 1:p-2
            for i = 1:p-1-j
              si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
              % ------------
              % origin @ vertex 1
              % ------------
              B(offset+1,:,1) = L3.*L4.*kernel3m2(:,i).*kernel4m2(:,j) + ...
                      L234.*(-dKernel3m2(:,i).*kernel4m2(:,j) - kernel3m2(:,i).*dKernel4m2(:,j));
              B(offset+1,:,2) = L2.*L4.*kernel3m2(:,i).*kernel4m2(:,j) + ...
                      L234.*(dKernel3m2(:,i).*kernel4m2(:,j) + 0*kernel3m2(:,i).*dKernel4m2(:,j));
              B(offset+1,:,3) = L2.*L3.*kernel3m2(:,i).*kernel4m2(:,j) + ...
                      L234.*(0*dKernel3m2(:,i).*kernel4m2(:,j) + kernel3m2(:,i).*dKernel4m2(:,j));
              % ------------
              % origin @ vertex 2
              % ------------          
              B(offset+2,:,1) = sj*(L3.*L4.*kernel4m3(:,i).*kernel3m2(:,j) + ...
                      L234.*(0*dKernel4m3(:,i).*kernel3m2(:,j) - kernel4m3(:,i).*dKernel3m2(:,j)));
              B(offset+2,:,2) = sj*(L2.*L4.*kernel4m3(:,i).*kernel3m2(:,j) + ...
                      L234.*(-dKernel4m3(:,i).*kernel3m2(:,j) + kernel4m3(:,i).*dKernel3m2(:,j)));
              B(offset+2,:,3) = sj*(L2.*L3.*kernel4m3(:,i).*kernel3m2(:,j) + ...
                      L234.*(dKernel4m3(:,i).*kernel3m2(:,j) + 0*kernel4m3(:,i).*dKernel3m2(:,j)));
              % ------------
              % origin @ vertex 3
              % ------------      
              B(offset+3,:,1) = sij*(L3.*L4.*kernel4m2(:,i).*kernel4m3(:,j) + ...
                      L234.*(-dKernel4m2(:,i).*kernel4m3(:,j) + 0*kernel4m2(:,i).*dKernel4m3(:,j)));
              B(offset+3,:,2) = sij*(L2.*L4.*kernel4m2(:,i).*kernel4m3(:,j) + ...
                      L234.*(0*dKernel4m2(:,i).*kernel4m3(:,j) - kernel4m2(:,i).*dKernel4m3(:,j)));
              B(offset+3,:,3) = sij*(L2.*L3.*kernel4m2(:,i).*kernel4m3(:,j) + ...
                      L234.*(dKernel4m2(:,i).*kernel4m3(:,j) + kernel4m2(:,i).*dKernel4m3(:,j)));
              offset = offset + 3;
            end
          end
          % face 4
          for j = 1:p-2
            for i = 1:p-1-j
              si = (-1)^(i-1); sj = (-1)^(j-1); sij = si*sj;
              % ------------
              % origin @ vertex 1
              % ------------
              B(offset+1,:,1) = -L3.*L4.*kernel3m1(:,i).*kernel4m1(:,j) + ...
                      L134.*(dKernel3m1(:,i).*kernel4m1(:,j) + kernel3m1(:,i).*dKernel4m1(:,j));
              B(offset+1,:,2) = L4.*(L1-L3).*kernel3m1(:,i).*kernel4m1(:,j) + ...
                      L134.*(2*dKernel3m1(:,i).*kernel4m1(:,j) + kernel3m1(:,i).*dKernel4m1(:,j));
              B(offset+1,:,3) = L3.*(L1-L4).*kernel3m1(:,i).*kernel4m1(:,j) + ...
                      L134.*(dKernel3m1(:,i).*kernel4m1(:,j) + 2*kernel3m1(:,i).*dKernel4m1(:,j));
              % ------------
              % origin @ vertex 2
              % ------------          
              B(offset+2,:,1) = sj*(-L3.*L4.*kernel4m3(:,i).*kernel3m1(:,j) + ...
                      L134.*(0*dKernel4m3(:,i).*kernel3m1(:,j) + kernel4m3(:,i).*dKernel3m1(:,j)));
              B(offset+2,:,2) = sj*(L4.*(L1-L3).*kernel4m3(:,i).*kernel3m1(:,j) + ...
                      L134.*(-dKernel4m3(:,i).*kernel3m1(:,j) + 2*kernel4m3(:,i).*dKernel3m1(:,j)));
              B(offset+2,:,3) = sj*(L3.*(L1-L4).*kernel4m3(:,i).*kernel3m1(:,j) + ...
                      L134.*(dKernel4m3(:,i).*kernel3m1(:,j) + kernel4m3(:,i).*dKernel3m1(:,j)));
              % ------------
              % origin @ vertex 3
              % ------------      
              B(offset+3,:,1) = sij*(-L3.*L4.*kernel4m1(:,i).*kernel4m3(:,j) + ...
                      L134.*(dKernel4m1(:,i).*kernel4m3(:,j) + 0*kernel4m1(:,i).*dKernel4m3(:,j)));
              B(offset+3,:,2) = sij*(L4.*(L1-L3).*kernel4m1(:,i).*kernel4m3(:,j) + ...
                      L134.*(dKernel4m1(:,i).*kernel4m3(:,j) - kernel4m1(:,i).*dKernel4m3(:,j)));
              B(offset+3,:,3) = sij*(L3.*(L1-L4).*kernel4m1(:,i).*kernel4m3(:,j) + ...
                      L134.*(2*dKernel4m1(:,i).*kernel4m3(:,j) + kernel4m1(:,i).*dKernel4m3(:,j)));
              offset = offset + 3;
            end
          end
          % inner
          L1234 = L1.*L234;
          for k = 1:p-3
            for j = 1:p-2-k
              for i = 1:p-1-j-k
                KKK = kernel2m1(:,i).*kernel3m1(:,j).*kernel4m1(:,k);
                S = dKernel2m1(:,i).*kernel3m1(:,j).*kernel4m1(:,k) + ...
                    kernel2m1(:,i).*dKernel3m1(:,j).*kernel4m1(:,k) + ...
                    kernel2m1(:,i).*kernel3m1(:,j).*dKernel4m1(:,k);
                B(offset+1,:,1) = L3.*L4.*(L1-L2).*KKK + L1234.* ...
                   (S + dKernel2m1(:,i).*kernel3m1(:,j).*kernel4m1(:,k));
                B(offset+1,:,2) = L2.*L4.*(L1-L3).*KKK + L1234.* ...
                   (S + kernel2m1(:,i).*dKernel3m1(:,j).*kernel4m1(:,k));
                B(offset+1,:,3) = L2.*L3.*(L1-L4).*KKK + L1234.* ...
                   (S + kernel2m1(:,i).*kernel3m1(:,j).*dKernel4m1(:,k));
                offset = offset + 1;
             end
            end
          end
      end
      B = permute(B, [1 2 4 3]);
    end
  end
  methods % kernel functions
    function [N,dN,d2N,d3N] = getKernel(obj, x)
      N = zeros(numel(x), obj.order-1);
      if nargout>1
        dN = zeros(numel(x),obj.order-1);
      end
      if nargout>2
        d2N = zeros(numel(x),obj.order-1);
      end
      if nargout>3
        d3N = zeros(numel(x),obj.order-1);
      end
      coeff = cell(20,1);
      coeff{1} = -2*sqrt(3/2);
      coeff{2} = -2*sqrt(5/2)*[1 0];
      coeff{3} = -0.5*sqrt(7/2)*[5 0 -1];
      coeff{4} = -0.5*sqrt(9/2)*[7 0 -3 0];
      coeff{5} = -0.25*sqrt(11/2)*[21 0 -14 0 1];
      coeff{6} = -0.25*sqrt(13/2)*[33 0 -30 0 5 0];
      coeff{7} = -1/32*sqrt(15/2)*[429 0 -495 0 135 0 -5];
      coeff{8} = -1/32*sqrt(17/2)*[715 0 -1001 0 385 0 -35 0];
      coeff{9} = -1/64*sqrt(19/2)*[2431 0 -4004 0 2002 0 -308 0 7];
      coeff{10} = -1/64*sqrt(21/2)*[4199 0 -7956 0 4914 0 -1092 0 63 0];
      coeff{11} = -1/256*sqrt(23/2)*[29393 0 -62985 0 46410 0 -13650 0 1365 0 -21];
      coeff{12} = -1/256*sqrt(25/2)*[52003 0 -124355 0 106590 0 -39270 0 5775 0 -231 0];
      coeff{13} = -1/512*sqrt(27/2)*[185725 0 -490314 0 479655 0 -213180 0 42075 0 -2970 0 33];
      coeff{14} = -1/512*sqrt(29/2)*[334305 0 -965770 0 1062347 0 -554268 0 138567 0 -14586 0 429 0];            
      coeff{15} = -1/8192*sqrt(31/2)*[9694845 0 -30421755 0 37182145 0 -22309287 0 6789783 0 -969969 0 51051 0 -429];
      coeff{16} = -1/8192*sqrt(33/2)*[17678835 0 -59879925 0 80528175 0 -54679625 0 19684665 0 -3594591 0 285285 0 -6435 0];
      coeff{17} = -1/16384*sqrt(35/2)*[64822395 0 -235717800 0 345972900 0 -262462200 0 109359250 0 -24496472 0 2662660 0 -108680 0 715];
      coeff{18} = -1/16384*sqrt(37/2)*[119409675 0 -463991880 0 738168900 0 -619109400 0 293543250 0 -78278200 0 10958948 0 -680680 0 12155 0];
      coeff{19} = -1/65536*sqrt(39/2)*[883631595 0 -3653936055 0 6263890380 0 -5757717420 0 3064591530 0 -951080130 0 164384220 0 -14090076 0 459459 0 -2431];
      coeff{20} = -1/65536*sqrt(41/2)*[1641030105 0 -7195285845 0 13223768580 0 -13223768580 0 7814045070 0 -2772725670 0 573667380 0 -63740820 0 3187041 0 -46189 0];
      for n = 0:obj.order-2
        %
        N(:,n+1) = polyval(coeff{n+1},x);
        if nargout>1
          dCoeff = coeff{n+1}(1:end-1).*(n:-1:1);
          dN(:,n+1) = polyval(dCoeff,x);
        end
        if nargout>2
          d2Coeff = dCoeff(1:end-1).*(n-1:-1:1);
          d2N(:,n+1) = polyval(d2Coeff,x);
        end
        if nargout>3
          d3Coeff = d2Coeff(1:end-1).*(n-2:-1:1);
          d3N(:,n+1) = polyval(d3Coeff,x);
        end
      end
    end
  end
end