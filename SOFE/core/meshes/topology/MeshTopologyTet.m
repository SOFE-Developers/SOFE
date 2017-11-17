classdef MeshTopologyTet < MeshTopology
  methods % constructor
    function obj = MeshTopologyTet(elem)
      obj = obj@MeshTopology(3);
      obj.update(elem);
      obj.isSimplex = 1;
    end
    function update(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{3,1} = [elem(:,[1,2,3]); elem(:,[1,2,4]); ...
                               elem(:,[2,3,4]); elem(:,[1 3 4])];
      [obj.connectivity{3,1}, ~, e2F] = unique(sort(obj.connectivity{3,1},2),'rows');    
      obj.connectivity{2,1} = [elem(:,[1,2]); elem(:,[2,3]); elem(:,[1,3]); ...
                               elem(:,[1 4]); elem(:,[2,4]); elem(:,[3 4])];
      [obj.connectivity{2,1}, ~, e2Ed] = unique(sort(obj.connectivity{2,1},2),'rows');    
      obj.connectivity{4,3} = reshape(e2F,[], 4);
      obj.connectivity{4,2} = reshape(e2Ed,[], 6);
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
      R = zeros(obj.getNumber(2), 3);
      e2F = obj.getElem2Face();
      e2E = obj.getElem2Edge();
      orientF = obj.getOrientation(3,2);
      [~, ind] = unique(e2F);
      [elems, type] = ind2sub([obj.getNumber(3), 4], ind);
      nodeIxAtFace = [1 2 3; 1 5 4; 2 6 5; 3 6 4];
      for t = 1:4
        for k = 0:2
          I = (type == t) & (abs(orientF(ind)) == k+1);
          R(I,:) = e2E(elems(I), circshift(nodeIxAtFace(t,:)',-k));
        end
      end
      I = orientF(ind)<0;
      R(I,:) = R(I, [3 2 1]);
    end
    function R = getOrientation(obj, dim1, dim2, varargin)
      R = [];
      switch dim2
        case 1
          switch dim1
            case 3
              e = obj.getEntity(3);
              R = ones(size(e,1), 6);
              R(e(:,1)>e(:,2),1) = -1;
              R(e(:,2)>e(:,3),2) = -1;
              R(e(:,1)>e(:,3),3) = -1;
              R(e(:,1)>e(:,4),4) = -1;
              R(e(:,2)>e(:,4),5) = -1;
              R(e(:,3)>e(:,4),6) = -1;
            case 2
              R = ones(obj.getNumber(2), 3);
            otherwise
              return
          end
        case 2
          switch dim1
            case 3
              e = obj.getEntity(3);
              face = [1 2 3; 1 2 4; 2 3 4; 1 3 4];
              R = zeros(size(e,1),4);
              for i = 1:4
                  [~, R(:,i)] = min(e(:,face(i,:)),[],2);
                  [~, P] = sort(e(:,face(i,:)),2);
                  even = (P(:,1)<P(:,2) & P(:,2)<P(:,3)) | ...
                         (P(:,2)<P(:,3) & P(:,3)<P(:,1)) | ...
                         (P(:,3)<P(:,1) & P(:,1)<P(:,2));
                  R(~even, i) = -R(~even, i);
              end
            otherwise
              return
          end
        otherwise
          return
      end
      if nargin > 3
        R = R(varargin{1},:,:);
      end
    end
    function R = getNormalOrientation(obj)
      R = sign(obj.getOrientation(3, 2)); % nEx4
      R(:,[1 4]) = -R(:,[1 4]);
    end
  end
  methods % refinement
    function P = uniformRefine(obj)
      edges = obj.getEntity(1);
      el = obj.getEntity(3);
      nEd = obj.getNumber(1); nN = obj.getNumber(0);
      P = [speye(nN); sparse(repmat((1:nEd)',1,2), edges, 0.5)];
      newIndices = (nN+1 : nN+obj.getNumber(1));
      el = [el newIndices(obj.connectivity{4,2})];
      el = [el(:,[1 5 7 8]); el(:,[2 6 5 9]); ...
            el(:,[3 7 6 10]); el(:,[4 9 8 10]); ...
            el(:,[5 6 7 9]); el(:,[5 7 8 9]); ...
            el(:,[7 8 9 10]); el(:,[7 9 6 10])];
      obj.update(el);
    end
  end
  methods(Static = true)
    function R = getQuadRule(order)
      R{4} = GaussPoint();
      R{3} = GaussInt(order);
      R{2} = GaussTri(order);
      R{1} = GaussTet(order);
    end
    function R = isFeasible(points)
      tol = 1e-12;
      R = (all(points>-tol, 2) & 1-sum(points,2)>-tol);
    end
    function R = getCenterLoc()
      R = [1 1 1]/4;
    end
    function R = renumber(node, elem)
      v1 = node(elem(:,2),:) - node(elem(:,1),:);
      v2 = node(elem(:,3),:) - node(elem(:,1),:);
      v3 = node(elem(:,4),:) - node(elem(:,1),:);
      I = ((v1(:,1).*v2(:,2).*v3(:,3) + v1(:,2).*v2(:,3).*v3(:,1)+v1(:,3).*v2(:,1).*v3(:,2)) ...
      - (v1(:,3).*v2(:,2).*v3(:,1)+v1(:,2).*v2(:,1).*v3(:,3)+v1(:,1).*v2(:,3).*v3(:,2)))/6;
      I = I<0;
      if any(I)
        fprintf('Elements renumbered!\n');
        elem(I,:) = elem(I, [2 1 3 4]);
      end
      R = elem;
    end
  end
end