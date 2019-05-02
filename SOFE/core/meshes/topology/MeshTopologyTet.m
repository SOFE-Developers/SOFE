classdef MeshTopologyTet < MeshTopology
  methods % constructor
    function obj = MeshTopologyTet(elem)
      obj = obj@MeshTopology(3);
      obj.update(elem);
      obj.isSimplex = 1;
      obj.nESub = [4 6 4 1];
      obj.nO = 6;
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
        for k = 1:3
          I = reshape(type == t,[],1) & (ceil(0.5*reshape(orientF(ind),[],1)) == k);
          R(I,:) = e2E(elems(I), circshift(nodeIxAtFace(t,:)',1-k));
        end
      end
      I = mod(orientF(ind),2)==0;
      R(I,:) = R(I, [3 2 1]);
    end
    function R = getOrientation(obj, dim, d, varargin) % [I]
      R = [];
      switch d
        case 1
          switch dim
            case 3
              e = obj.getEntity(3, varargin{:});
              R = ones(size(e,1), 6);
              R(e(:,1)>e(:,2),1) = 2;
              R(e(:,2)>e(:,3),2) = 2;
              R(e(:,1)>e(:,3),3) = 2;
              R(e(:,1)>e(:,4),4) = 2;
              R(e(:,2)>e(:,4),5) = 2;
              R(e(:,3)>e(:,4),6) = 2;
            otherwise
              % always positive
              return
          end
        case 2
          switch dim
            case 3
              e = obj.getEntity(3, varargin{:});
              face = [1 2 3; 1 2 4; 2 3 4; 1 3 4];
              R = zeros(size(e,1),4);
              for i = 1:4
                [~, R(:,i)] = min(e(:,face(i,:)),[],2);
                [~, P] = sort(e(:,face(i,:)),2);
                even = (P(:,1)<P(:,2) & P(:,2)<P(:,3)) | ...
                       (P(:,2)<P(:,3) & P(:,3)<P(:,1)) | ...
                       (P(:,3)<P(:,1) & P(:,1)<P(:,2));
                R(:,i) = 2*R(:,i) - 1;
                R(~even, i) = R(~even, i) + 1;
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
      R(:,[1 4]) = -R(:,[1 4]);
    end
  end
  methods % refinement
    function R = uniformRefine(obj)
      ed = obj.getEntity(1);
      el = obj.getEntity(3);
      nN = obj.getNumber(0); nEd = obj.getNumber(1);
      edRange = (1:nEd)';
      %
      R = [speye(nN); fsparse(repmat((1:nEd)',1,2), ed, 0.5)];
      %
      newIndices = nN + edRange;
      el = [el newIndices(obj.connectivity{4,2})];
      el = [el(:,[1 5 7 8]); el(:,[5 2 6 9]); ...
            el(:,[7 6 3 10]); el(:,[8 9 10 4]); ...
            el(:,[6 7 5 9]); el(:,[5 7 8 9]); ...
            el(:,[8 10 9 7]); el(:,[7 9 6 10])];
            keyboard
      obj.update(el);
    end
    function R = uniformRefineFast(obj)
      ed = obj.getEntity(1);
      fc = obj.getEntity(2);
      el = obj.getEntity(3);
      e2Ed = obj.connectivity{4,2};
      f2Ed = obj.connectivity{3,2};
%      oF = obj.getOrientation(3,2)==???;
      nN = obj.getNumber(0); nEd = obj.getNumber(1);
      nF = obj.getNumber(2); nE = obj.getNumber(3);
      edRange = (1:nEd)';
      fcRange = (1:nF)';
      elRange = (1:nE)';
      %
      R = [speye(nN); fsparse(repmat((1:nEd)',1,2), ed, 0.5)];
      %
      keyboard
      newIndices = nN + edRange;
      el = [el newIndices(e2Ed)];
      el = [el(:,[1 5 7 8]); el(:,[5 2 6 9]); ...
            el(:,[7 6 3 10]); el(:,[8 9 10 4]); ...
            el(:,[6 7 5 9]); el(:,[5 7 8 9]); ...
            el(:,[8 10 9 7]); el(:,[7 9 6 10])];
      fc = [fc newIndices(f2Ed)];
      fc = [[fc(:,1) sort(fc(:,[4 6]),2)]; ...
            [fc(:,2) sort(fc(:,[4 5]),2)]; ...
            [fc(:,3) sort(fc(:,[5 6]),2)]; ...
            sort(fc(:,[4 5 6]),2)];
            % copy & adapt of tri, TODO: check + 4 inner faces
      ed = [ed(:,1) newIndices; ed(:,2) newIndices; ...
            sort([nN+f2Ed(:,[1 2]);nN+f2Ed(:,[2 3]);nN+f2Ed(:,[1 3]);nN+e2Ed(:,[3 5])],2)];
      %
      f2Ed = [[f2Ed(:,1), 2*(nEd+nF)+fcRange, f2Ed(:,3)]; ...
              [nEd+f2Ed(:,1), f2Ed(:,2), 2*nEd+fcRange]; ...
              [(2*nEd+nF)+fcRange, nEd+f2Ed(:,2), nEd+f2Ed(:,3)]; ...
              [(2*nEd+nF)+fcRange, 2*(nEd+nF)+fcRange, 2*nEd+fcRange]];;
              % copy & adapt of tri, TODO: check + 4 inner faces
      e2F = []; % TODO
      e2Ed = []; % TODO
      %
      obj.connectivity{2,1} = ed;
      obj.connectivity{3,1} = fc;
      obj.connectivity{4,1} = el;
      obj.connectivity{3,2} = f2Ed;
      obj.connectivity{4,3} = e2F;
      obj.connectivity{4,2} = e2Ed;
      obj.connectivity{1,1} = (1:size(R,1))';
      obj.connectivity{2,2} = (1:size(ed,1))';
      obj.connectivity{3,3} = (1:size(fc,1))';
      obj.connectivity{4,4} = (1:size(el,1))';
    end
  end
  methods(Static = true)
    function R = isFeasible(points, varargin) % [tol
      if ~isempty(varargin), tol = varargin{1}; else, tol = 1e-12; end
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
    function R = upliftPoints(points, fLoc, orient)
      switch orient
        case 1
          B = [1 0;0 1]; b = [0;0];
        case 2
          B = [0 1;1 0]; b = [0;0];
        case 3
          B = [-1 -1;1 0]; b = [1;0];
        case 4
          B = [-1 -1;0 1]; b = [1;0];
        case 5
          B = [0 1;-1 -1]; b = [0;1];
        case 6
          B = [1 0;-1 -1]; b = [0;1];
      end
      R = bsxfun(@plus, B*points', b)';
      switch fLoc
        case 1
          B = [1 0;0 1;0 0]; b = [0;0;0];
        case 2
          B = [1 0;0 0;0 1]; b = [0;0;0];
        case 3
          B = [-1 -1;1 0;0 1]; b = [1;0;0];
        case 4
          B = [0 0;1 0;0 1]; b = [0;0;0];
      end
      R = bsxfun(@plus, B*R', b)';
    end
  end
end