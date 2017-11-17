classdef MeshTopologyTri < MeshTopology
  methods % constructor
    function obj = MeshTopologyTri(elem)
      obj = obj@MeshTopology(2);
      obj.update(elem);
      obj.isSimplex = 1;
    end
    function update(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{2,1} = [elem(:,[1,2]); elem(:,[2,3]); elem(:,[1,3])];
      [obj.connectivity{2,1}, ~, e2F] = unique(sort(obj.connectivity{2,1},2),'rows');      
      obj.connectivity{3,2} = reshape(e2F, size(elem,1), []);
      %
      obj.connectivity{1,1} = (1:max(obj.connectivity{2,1}(:)))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
    end
  end
  methods % connectivity information
    function R = getOrientation(obj, varargin)
      e = obj.getEntity('0');
      R = ones(size(e));
      R(e(:,1)>e(:,2),1) = -1;
      R(e(:,2)>e(:,3),2) = -1;
      R(e(:,1)>e(:,3),3) = -1;
      if nargin > 3
        R = R(varargin{3},:);
      end
    end
    function R = getNormalOrientation(obj, varargin)
      R = obj.getOrientation();
      R(:,3) = -R(:,3);
    end
  end
  methods % refinement & manipulation
    function P = uniformRefine(obj)
      el = obj.getEntity(2);
      nF = obj.getNumber(1); nN = obj.getNumber(0);
      faces = obj.getEntity(1);
      P = [speye(nN); sparse(repmat((1:nF)',1,2), faces, 0.5)];
      newIndices = nN + (1:nF);
      el = [el newIndices(obj.connectivity{3,2})];
      el = [el(:,[1 4 6]);el(:,[4 2 5]);el(:,[6 5 3]);el(:,[5 6 4])];
      obj.update(el);
    end
    function flipFace(obj, I)
      [f2e, type] = obj.getFace2Elem();
      elem = obj.getEntity(2);
      face = obj.getEntity(1,I);
      idxL = elem(f2e(I,1) + size(elem,1)*mod(type(I,1)-1+2,3));
      idxR = elem(f2e(I,2) + size(elem,1)*mod(type(I,2)-1+2,3));
      elem(f2e(I,1),:) = [idxR, idxL, face(:,1)];
      elem(f2e(I,2),:) = [idxL, idxR, face(:,2)];
      obj.update(elem);
    end
  end
  methods(Static = true)
    function R = renumber(node, elem)
      % positive jacobian
      if isempty(elem), R = []; return; end
      v1 = node(elem(:,2),:) - node(elem(:,1),:);
      v2 = node(elem(:,3),:) - node(elem(:,1),:);
      I = v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1);
      I = I<0;
      if any(I)
        fprintf('Elements renumbered!\n');
        elem(I,:) = elem(I, [2 1 3]);
      end
      R = elem;
    end
    function R = isFeasible(points)
      tol = 1e-12;
      R = (all(points>-tol, 2) & 1-sum(points,2)>-tol);
    end
    function R = getCenterLoc()
      R = [1 1]/3;
    end
    function R = getQuadRule(quadOrder)
      R{3} = GaussPoint();
      R{2} = GaussInt(quadOrder);
      R{1} = GaussTri(quadOrder);
    end
    function R = upliftPoints(points, fLoc, orient)
      zz = zeros(size(points));
      if orient<0
        points = 1-points;
      end
      switch fLoc
        case 1
          R = [points, zz];
        case 2
          R = [1-points, points];
        case 3
          R = [zz, points];
      end
    end
  end
end