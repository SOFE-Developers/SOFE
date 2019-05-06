classdef MeshTopologyHex < MeshTopology
  methods % constructor
    function obj = MeshTopologyHex(elem)
      obj = obj@MeshTopology(3);
      obj.update(elem);
      obj.isSimplex = 0;
      obj.nESub = [8 12 6 1];
      obj.nO = 8;
    end
    function update(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{3,1} = [elem(:,[1,2,3,4]); elem(:,[5,6,7,8]); ...
                               elem(:,[1,3,5,7]); elem(:,[2,4,6,8]); ...
                               elem(:,[1,2,5,6]); elem(:,[3,4,7,8])];
      [~,I] = min(obj.connectivity{3,1},[],2);
      obj.connectivity{3,1}(I==1,:) = obj.connectivity{3,1}(I==1,[1 2 3 4]);
      obj.connectivity{3,1}(I==2,:) = obj.connectivity{3,1}(I==2,[2 1 4 3]);
      obj.connectivity{3,1}(I==3,:) = obj.connectivity{3,1}(I==3,[3 1 4 2]);
      obj.connectivity{3,1}(I==4,:) = obj.connectivity{3,1}(I==4,[4 2 3 1]);
      obj.connectivity{3,1}(:,2:3) = sort(obj.connectivity{3,1}(:,2:3),2);
      [obj.connectivity{3,1}, ~, e2F] = unique(obj.connectivity{3,1}, 'rows');    
      obj.connectivity{2,1} = [elem(:,[1,2]); elem(:,[3,4]); elem(:,[5,6]); elem(:,[7,8]); ...
                               elem(:,[1,3]); elem(:,[5,7]); elem(:,[2,4]); elem(:,[6,8]); ...
                               elem(:,[1,5]); elem(:,[2,6]); elem(:,[3,7]); elem(:,[4,8])];
      [obj.connectivity{2,1}, ~, e2Ed] = unique(sort(obj.connectivity{2,1},2),'rows');    
      obj.connectivity{4,3} = reshape(e2F, [], 6);
      obj.connectivity{4,2} = reshape(e2Ed, [], 12);
      %
      obj.connectivity{1,1} = (1:max(obj.connectivity{3,1}(:)))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
      obj.connectivity{4,4} = (1:size(obj.connectivity{4,1},1))';
      %
      obj.connectivity{3,2} = obj.getFace2Edge();
    end
  end
  methods % connectivity information
    function R = getElem2Edge(obj)
      R = obj.connectivity{4,2};
    end
    function R = getFace2Edge(obj)
      R = zeros(obj.getNumber(2), 4);
      e2F = obj.getElem2Face();
      e2E = obj.getElem2Edge();
      orientF = obj.getOrientation(3,2);
      [~, ind] = unique(e2F);
      orientF = reshape(orientF(ind), size(ind));
      [father, type] = ind2sub([obj.getNumber(3), 6], ind);
      nodeIxAtFace = [1 2 5 7; 3 4 6 8; 5 6 9 11; 7 8 10 12; 1 3 9 10; 2 4 11 12];
      flags = [1 1 1;1 1 -1; 1 -1 -1; -1 1 1; -1 1 -1; 1 -1 1; -1 -1 1; -1 -1 -1];
      for t = 1:6
        for k = 1:8 % to check, compare with getDoFEnum@LagrangeElement
          I = (type == t) & (orientF == k);
          idxMap = [1 3; 2 4];
          if flags(k,3)<0, idxMap = idxMap(:,[2 1]); end
          if flags(k,1)<0, idxMap(:,2) = idxMap([2 1],2); end
          if flags(k,2)<0, idxMap(:,1) = idxMap([2 1],1); end
          R(I,:) = e2E(father(I), nodeIxAtFace(t,idxMap(:)));
        end
      end
    end
    function R = getOrientation(obj, dim, d, varargin) % [I]
      R = [];
      switch d
        case 1
          switch dim
            case 3
              e = obj.getEntity(3, varargin{:});   
              R = ones(size(e,1), 12);
              R(e(:,1)>e(:,2),1) = 2;
              R(e(:,3)>e(:,4),2) = 2;
              R(e(:,5)>e(:,6),3) = 2;
              R(e(:,7)>e(:,8),4) = 2;
              R(e(:,1)>e(:,3),5) = 2;
              R(e(:,5)>e(:,7),6) = 2;
              R(e(:,2)>e(:,4),7) = 2;
              R(e(:,6)>e(:,8),8) = 2;
              R(e(:,1)>e(:,5),9) = 2;
              R(e(:,2)>e(:,6),10) = 2;
              R(e(:,3)>e(:,7),11) = 2;
              R(e(:,4)>e(:,8),12) = 2;
            case 2
              e = obj.getEntity(2, varargin{:});  
              R = ones(size(e));
              R(e(:,1)>e(:,2),1) = 2;
              R(e(:,3)>e(:,4),2) = 2;
              R(e(:,1)>e(:,3),3) = 2;
              R(e(:,2)>e(:,4),4) = 2;
            otherwise
              return
          end
        case 2
          switch dim
            case 3
              e = obj.getEntity(3, varargin{:});
              R = ones(size(e,1), 6);
              face = [1 2 3 4;
                      5 6 7 8;
                      1 3 5 7;
                      2 4 6 8;
                      1 2 5 6;
                      3 4 7 8];
              for i = 1:6
                [~, minVx] = min(e(:,face(i,:)),[],2);
                R((minVx == 1) & (e(:,face(i,2)) < e(:,face(i,3))),i) = 1;
                R((minVx == 1) & (e(:,face(i,2)) > e(:,face(i,3))),i) = 2;
                R((minVx == 2) & (e(:,face(i,1)) > e(:,face(i,4))),i) = 3;
                R((minVx == 2) & (e(:,face(i,1)) < e(:,face(i,4))),i) = 4;
                R((minVx == 3) & (e(:,face(i,1)) < e(:,face(i,4))),i) = 5;
                R((minVx == 3) & (e(:,face(i,1)) > e(:,face(i,4))),i) = 6;
                R((minVx == 4) & (e(:,face(i,2)) > e(:,face(i,3))),i) = 7;
                R((minVx == 4) & (e(:,face(i,2)) < e(:,face(i,3))),i) = 8;
              end
            otherwise
              return
          end
        otherwise
          return
      end
    end
    function R = getNormalOrientation(obj, varargin) % [I]
      R = 2*mod(obj.getOrientation(3, 2, varargin{:}),2)-1; % nExnF
      R(:,[1 3 6]) = -R(:,[1 3 6]);
    end
  end
  methods % refine
    function R = uniformRefine(obj)
      el = obj.getEntity(3);
      faces = obj.getEntity(2);
      edges = obj.getEntity(1);
      nN = obj.getNumber(0); nEd = obj.getNumber(1);
      nF = obj.getNumber(2); nE = obj.getNumber(3);
      R = [speye(nN); fsparse(repmat((1:nEd)',1,2), edges, 0.5); ...
                      fsparse(repmat((1:nF)',1,4), faces, 0.25); ...
                      fsparse(repmat((1:nE)',1,8), el, 0.125)];
      newIndicesEd = (nN+1 : nN+nEd);
      newIndicesF = (nN+nEd+1 : nN+nEd+nF);
      newIndicesE = (nN+nEd+nF+1 : nN+nEd+nF+nE)';
      el = [el,newIndicesEd(obj.connectivity{4,2}),newIndicesF(obj.connectivity{4,3}),newIndicesE];
      el = [el(:,[1 9 13 21 17 25 23 27]); el(:,[9 2 21 15 25 18 27 24]); ...
            el(:,[13 21 3 10 23 27 19 26]); el(:,[21 15 10 4 27 24 26 20]); ...
            el(:,[17 25 23 27 5 11 14 22]); el(:,[25 18 27 24 11 6 22 16]); ...
            el(:,[23 27 19 26 14 22 7 12]); el(:,[27 24 26 20 22 16 12 8])];
      obj.update(el);
    end
    function R = uniformRefineFast(obj)
      R = obj.uniformRefine(); % TODO
    end
  end
  methods(Static = true)
    function R = isFeasible(points, varargin) % [tol
      if ~isempty(varargin), tol = varargin{1}; else, tol = 1e-12; end
      R = (points(:,1)>-tol & points(:,1)<1+tol & ...
           points(:,2)>-tol & points(:,2)<1+tol & ...
           points(:,3)>-tol & points(:,3)<1+tol );
    end
    function R = getCenterLoc()
      R = [1 1 1]/2;
    end
    function R = upliftPoints(points, fLoc, orient)
      switch orient
        case 1
          B = [1 0;0 1]; b = [0;0];
        case 2
          B = [0 1;1 0]; b = [0;0];
        case 3
          B = [0 -1;1 0]; b = [1;0];
        case 4
          B = [-1 0;0 1]; b = [1;0];
        case 5
          B = [0 1;-1 0]; b = [0;1];
        case 6
          B = [1 0;0 -1]; b = [0;1];
        case 7
          B = [-1 0;0 -1]; b = [1;1];
        case 8
          B = [0 -1;-1 0]; b = [1;1];
      end
      R = bsxfun(@plus, B*points', b)';
      switch fLoc
        case 1
          B = [1 0;0 1;0 0]; b = [0;0;0];
        case 2
          B = [1 0;0 1;0 0]; b = [0;0;1];
        case 3
          B = [0 0;1 0;0 1]; b = [0;0;0];
        case 4
          B = [0 0;1 0;0 1]; b = [1;0;0];
        case 5
          B = [1 0;0 0;0 1]; b = [0;0;0];
        case 6
          B = [1 0;0 0;0 1]; b = [0;1;0];
      end
      R = bsxfun(@plus, B*R', b)';
    end
  end
end