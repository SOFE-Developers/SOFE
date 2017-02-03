classdef MeshTopologyHex < MeshTopology
  methods % constructor
    function obj = MeshTopologyHex(nodes, elems, dimP)
      obj = obj@MeshTopology(nodes, elems, dimP);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{3,1} = [obj.connectivity{4,1}(:,[1,2,3,4]); ...
                               obj.connectivity{4,1}(:,[5,6,7,8]); ...
                               obj.connectivity{4,1}(:,[1,3,5,7]); ...
                               obj.connectivity{4,1}(:,[2,4,6,8]); ...
                               obj.connectivity{4,1}(:,[1,2,5,6]); ...
                               obj.connectivity{4,1}(:,[3,4,7,8])];
      [~,I] = min(obj.connectivity{3,1},[],2);
      obj.connectivity{3,1}(I==1,:) = obj.connectivity{3,1}(I==1,[1 2 3 4]);
      obj.connectivity{3,1}(I==2,:) = obj.connectivity{3,1}(I==2,[2 1 4 3]);
      obj.connectivity{3,1}(I==3,:) = obj.connectivity{3,1}(I==3,[3 1 4 2]);
      obj.connectivity{3,1}(I==4,:) = obj.connectivity{3,1}(I==4,[4 2 3 1]);
      obj.connectivity{3,1}(:,2:3) = sort(obj.connectivity{3,1}(:,2:3),2);
      [obj.connectivity{3,1}, ~, e2F] = unique(obj.connectivity{3,1}, 'rows');    
      obj.connectivity{2,1} = [obj.connectivity{4,1}(:,[1,2]); obj.connectivity{4,1}(:,[3,4]); ...
                               obj.connectivity{4,1}(:,[5,6]); obj.connectivity{4,1}(:,[7,8]); ...
                               obj.connectivity{4,1}(:,[1,3]); obj.connectivity{4,1}(:,[5,7]); ...
                               obj.connectivity{4,1}(:,[2,4]); obj.connectivity{4,1}(:,[6,8]); ...
                               obj.connectivity{4,1}(:,[1,5]); obj.connectivity{4,1}(:,[2,6]); ...
                               obj.connectivity{4,1}(:,[3,7]); obj.connectivity{4,1}(:,[4,8])];
      [obj.connectivity{2,1}, ~, e2Ed] = unique(sort(obj.connectivity{2,1},2),'rows');    
      obj.connectivity{4,3} = reshape(e2F, [], 6);
      obj.connectivity{4,2} = reshape(e2Ed, [], 12);
      obj.connectivity{3,2} = obj.getFace2Edge();
      %
      obj.connectivity{1,1} = (1:size(obj.nodes,1))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
      obj.connectivity{4,4} = (1:size(obj.connectivity{4,1},1))';
    end
  end
  methods % connectivity information
    function R = getElem2Face(obj)
      R = obj.connectivity{4,3};
    end
    function R = getElem2Edge(obj)
      R = obj.connectivity{4,2};
    end
    function R = getFace2Edge(obj)
      R = zeros(obj.getNumber(2), 4);
      e2F = obj.getElem2Face();
      e2E = obj.getElem2Edge();
      orientF = obj.getOrientation(3,2);
      orient1 = orientF(:,:,1);
      orient2 = orientF(:,:,2);
      orient3 = orientF(:,:,3);
      [~, ind] = unique(e2F);
      orient1 = orient1(ind); orient2 = orient2(ind); orient3 = orient3(ind);
      [father, type] = ind2sub([obj.getNumber(3), 6], ind);
      nodeIxAtFace = [1 2 5 7; 3 4 6 8; 5 6 9 11; 7 8 10 12; 1 3 9 10; 2 4 11 12];
      for t = 1:6
        for s1 = -1:2:1
          for s2 = -1:2:1
            for k = -1:2:1
              I = (type == t) & (orient1 == s1) & (orient2 == s2) & (orient3 == k);
              idxMap = [1 3; 2 4];
              if k<0, idxMap = idxMap(:,[2 1]); end
              if s1<0, idxMap(:,2) = idxMap([2 1],2); end
              if s2<0, idxMap(:,1) = idxMap([2 1],1); end
              R(I,:) = e2E(father(I), nodeIxAtFace(t,idxMap(:)));
            end
          end
        end
      end
    end
    function R = getOrientation(obj, dim1, dim2)
      R = [];
      switch dim2
        case 1
          switch dim1
            case 3
              e = obj.getEntity(3);   
              R = ones(size(e,1), 12);
              R(e(:,1)>e(:,2),1) = -1;
              R(e(:,3)>e(:,4),2) = -1;
              R(e(:,5)>e(:,6),3) = -1;
              R(e(:,7)>e(:,8),4) = -1;
              R(e(:,1)>e(:,3),5) = -1;
              R(e(:,5)>e(:,7),6) = -1;
              R(e(:,2)>e(:,4),7) = -1;
              R(e(:,6)>e(:,8),8) = -1;
              R(e(:,1)>e(:,5),9) = -1;
              R(e(:,2)>e(:,6),10) = -1;
              R(e(:,3)>e(:,7),11) = -1;
              R(e(:,4)>e(:,8),12) = -1;
            case 2
              e = obj.getEntity(2);  
              R = ones(size(e));
              R(e(:,1)>e(:,2),1) = -1;
              R(e(:,3)>e(:,4),2) = -1;
              R(e(:,1)>e(:,3),3) = -1;
              R(e(:,2)>e(:,4),4) = -1;
          end
        case 2
          R = ones(3, obj.getNumber(3), 6); % ! three orient flags
          e = obj.getEntity(3); 
          face = [1 2 3 4;
                  5 6 7 8;
                  1 3 5 7;
                  2 4 6 8;
                  1 2 5 6;
                  3 4 7 8];
          for i = 1:6
            [~, minVx] = min(e(:,face(i,:)),[],2);
            I1p = (minVx == 1) & (e(:,face(i,2)) < e(:,face(i,3)));
            I1n = (minVx == 1) & (e(:,face(i,2)) > e(:,face(i,3)));
            I2n = (minVx == 2) & (e(:,face(i,1)) < e(:,face(i,4)));
            I2p = (minVx == 2) & (e(:,face(i,1)) > e(:,face(i,4)));
            I3p = (minVx == 3) & (e(:,face(i,1)) < e(:,face(i,4)));
            I3n = (minVx == 3) & (e(:,face(i,1)) > e(:,face(i,4)));
            I4n = (minVx == 4) & (e(:,face(i,2)) < e(:,face(i,3)));
            I4p = (minVx == 4) & (e(:,face(i,2)) > e(:,face(i,3)));
            R(:,I1p,i) = repmat([1 1 1]',   1, sum(I1p));
            R(:,I1n,i) = repmat([1 1 -1]',  1, sum(I1n));
            R(:,I2n,i) = repmat([-1 1 1]',  1, sum(I2n));
            R(:,I2p,i) = repmat([1 -1 -1]', 1, sum(I2p));
            R(:,I3p,i) = repmat([-1 1 -1]', 1, sum(I3p));
            R(:,I3n,i) = repmat([1 -1 1]',  1, sum(I3n));
            R(:,I4n,i) = repmat([-1 -1 -1]',1, sum(I4n));
            R(:,I4p,i) = repmat([-1 -1 1]', 1, sum(I4p));
          end
          R = permute(R,[2 3 1]);
      end
    end
    function R = getNormalOrientation(obj, varargin)
      orient = prod(obj.getOrientation([], 2), 3); % nExnF
      R = ones(obj.getNumber(3), 6); % nExnF
      R(:,[2 4 6]) =  orient(:,[2 4 6]);
      R(:,[1 3 5]) = -orient(:,[1 3 5]);
    end
  end
  methods % mesh information
    function R = getQuadRule(obj, order)
      R{4} = GaussPoint();
      R{3} = GaussInt(order);
      R{2} = GaussQuad(order);
      R{1} = GaussHex(order);
    end
    function R = getBarycenterRef(obj)
      R = [1 1 1]/2;
    end
    function R = getMeasure(obj, dim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      ee = obj.getEntity(dim); ee = ee(I,:);
      switch dim
        case 3 % element
          
        case 2 % face
          v1 = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          v2 = obj.nodes(ee(:,3),:) - obj.nodes(ee(:,1),:);
          R = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1));
          if all(all(obj.nodes(ee(:,1),:) + v1+v2 - obj.nodes(ee(:,4),:)))
            warning('! Area only valid for parallelograms !');
          end
        case 1 % edge
          v = obj.nodes(ee(:,2),:) - obj.nodes(ee(:,1),:);
          R = sum(v.^2,2).^0.5;
       end
    end
    function R = isFeasible(obj, points)
      tol = 1e-12;
      R = (points(:,1)>-tol & points(:,1)<1+tol & ...
           points(:,2)>-tol & points(:,2)<1+tol & ...
           points(:,3)>-tol & points(:,3)<1+tol );
    end
  end
  methods % refine
    function uniformRefine(obj)
      el = obj.getEntity(3);
      nN = obj.getNumber(0); nEd = obj.getNumber(1);
      nF = obj.getNumber(2); nE = obj.getNumber(3);
      % node coords
      obj.nodes = [obj.nodes; obj.getCenter(2); obj.getCenter(1); obj.getCenter(0)];
      % new node indices
      newNodeIxEd = (nN+1 : nN+nEd);
      newNodeIxF = (nN+nEd+1 : nN+nEd+nF);
      newNodeIxE = (nN+nEd+nF+1 : nN+nEd+nF+nE)';
      el = [el,newNodeIxEd(obj.connectivity{4,2}),newNodeIxF(obj.connectivity{4,3}),newNodeIxE];
      % elem
      obj.connectivity{4,1} = [el(:,[1 9 13 21 17 25 23 27]); el(:,[9 2 21 15 25 18 27 24]); ...
                               el(:,[13 21 3 10 23 27 19 26]); el(:,[21 15 10 4 27 24 26 20]); ...
                               el(:,[17 25 23 27 5 11 14 22]); el(:,[25 18 27 24 11 6 22 16]); ...
                               el(:,[23 27 19 26 14 22 7 12]); el(:,[27 24 26 20 22 16 12 8])];
      % faces
      obj.updateConnectivity();
      obj.notifyObservers();
    end
  end
  methods % display
    function show(obj, varargin)
      c = caxis();
      fc = obj.getEntity(2);
      I = obj.isSurface(varargin{:});
      h = trimesh(fc(I,[1 2 4 3]), obj.nodes(:,1), obj.nodes(:,2), obj.nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
      axis equal, axis tight, caxis(c);
    end
    function showEntity(obj, dim)
      center = obj.getCenter(dim);
      if dim == 0
        nE = size(obj.nodes, 1);
      else
        nE = obj.getNumber(dim);
      end
      switch dim
        case 3
          color = [1 0 1];
        case 2
          color = [0.1 0.1 0.1];
        case 1
          color = [0 0.5 1];
        case 0
          color = [0 1 0.5];
      end
      text(center(:,1), center(:,2), center(:,3), num2str((1:nE)'),'Color',color,'FontSize', 18);
    end
    function showNodeVector(obj, U, varargin)
      fc = obj.getEntity(2);
      if nargin > 2
        I = obj.isSurface(varargin{1});
      else
        I = obj.isBoundary();
      end
      h = trimesh(fc(I,[1 2 4 3]),obj.nodes(:,1),obj.nodes(:,2),obj.nodes(:,3), U(1:obj.getNumber(0)));
      set(h,'facecolor','interp','edgecolor','k');
      axis equal, axis tight
    end
  end
end