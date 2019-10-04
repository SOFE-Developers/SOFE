classdef MeshTopology < SOFE
  properties
    dimP
    connectivity
    isSimplex
    nESub, nO
  end
  methods % constructor
    function obj = MeshTopology(dimP)
      obj.dimP = dimP;
    end
  end
  methods(Static = true)
    function R = create(nodes, elem, dimP)
      switch size(elem, 2)
        case 1
          R = MeshTopologyPoint(elem);
        case 2
          R = MeshTopologyInt(elem);
        case 3
          %elem = MeshTopologyTri.renumber(nodes, elem);
          R = MeshTopologyTri(elem);
        case 4
          if dimP == 2
           R = MeshTopologyQuad(elem);
          else
            elem = MeshTopologyTet.renumber(nodes, elem);
            R = MeshTopologyTet(elem);
          end
        case 8
          R = MeshTopologyHex(elem);
      end
    end
  end
  methods % mesh information
    function R = getEntity(obj, dim, varargin) % [I]
      I = ':'; if nargin > 2, I = varargin{1}; end
      if ischar(dim), dim  = obj.dimP - str2double(dim); end % dim to codim
      R = obj.connectivity{dim+1}(I,:);
    end
    function R = getNumber(obj, varargin) % [dim]
      if isempty(varargin)
        R = zeros(obj.dimP+1,1);
        for d = 1:numel(R)
          R(d) = size(obj.connectivity{d, 1},1);
        end
      else
        dim = varargin{1};
        if ischar(dim), dim  = obj.dimP - str2double(dim); end % dim to codim
        R = size(obj.connectivity{dim+1, 1}, 1);
      end
    end
    function R = isBoundary(obj)
      e2F = obj.getElem2Face(); % nExnF
      R = accumarray(reshape(e2F(e2F>0),[],1), 1, [obj.getNumber('1') 1])==1; % nFx1
    end
    function R = getBoundary(obj)
      R = obj.getEntity('1', obj.isBoundary());
    end
    function R = isBoundaryNode(obj)
      R = unique(obj.getBoundary());
      R = accumarray(R,1,[obj.getNumber(0) 1])>0;
    end
  end
  methods % connectivity information
    function R = getElem2Face(obj, varargin)
      R = obj.connectivity{obj.dimP+1,obj.dimP};
      if nargin > 2
        R = R(varargin{1},:);
      end
    end
    function [R, fType, fOrient] = getFace2Elem(obj)
      if isempty(obj.connectivity{obj.dimP,obj.dimP+1})        
        nE = obj.getNumber('0'); nF = obj.getNumber('1');
        orient = 0.5*(3-obj.getNormalOrientation());
        e2F = obj.getElem2Face();
        obj.connectivity{obj.dimP,obj.dimP+1}{1} = accumarray([e2F(:), orient(:)], repmat((1:nE)',size(orient,2),1),[nF 2]);
        obj.connectivity{obj.dimP,obj.dimP+1}{2} = accumarray([e2F(:), orient(:)], kron((1:size(orient,2))',ones(nE,1)),[nF 2]);
        obj.connectivity{obj.dimP,obj.dimP+1}{3} = zeros(size(obj.connectivity{obj.dimP,obj.dimP+1}{1}));
        orientE = obj.getOrientation(obj.dimP, obj.dimP-1);
        for i = 1:2
          I = obj.connectivity{obj.dimP,obj.dimP+1}{1}(:,i) > 0;
          idx = obj.connectivity{obj.dimP,obj.dimP+1}{1}(I,i) + ...
                size(orientE,1)*(obj.connectivity{obj.dimP,obj.dimP+1}{2}(I,i)-1);
          obj.connectivity{obj.dimP,obj.dimP+1}{3}(I,i) = orientE(idx);
        end
      end
      R = obj.connectivity{obj.dimP,obj.dimP+1}{1};
      if nargout>1
        fType = obj.connectivity{obj.dimP,obj.dimP+1}{2};
      end
      if nargout>2
        fOrient = obj.connectivity{obj.dimP,obj.dimP+1}{3};
      end
    end
    function [R, nType] = getConnect(obj, dimFrom, dimTo)
      if dimFrom < dimTo
        if isempty(obj.connectivity{dimFrom+1,dimTo+1}) || nargout>1
          entity = obj.connectivity{dimTo+1, dimFrom+1};
          sz = size(entity);
          [entity,I] = sort(entity(:));
          count = accumarray(entity(:),1);
          maxCount = max(count);
          upperTri = triu(repmat((1:maxCount)',1,maxCount));
          count = upperTri(:,count);
          count(count==0) = [];
          idx = repmat((1:sz(1))', 1, sz(2));
          obj.connectivity{dimFrom+1,dimTo+1} = accumarray([entity, count(:)], idx(I));
          if nargout>1
            idx = kron(ones(sz(1),1), (1:sz(2)));
            nType = accumarray([entity, count(:)], idx(I));
          end
        end
      end
      R = obj.connectivity{dimFrom+1,dimTo+1};
    end
    function R = getConnect2(obj, dimFrom, dimTo)
      % slow version (3x) by sparse and transpose (just for validation)
      if dimFrom >= dimTo
        R = obj.connectivity{dimFrom+1, dimTo+1};
      else
        if isempty(obj.connectivity{dimFrom+1,dimTo+1})
          entity = obj.connectivity{dimTo+1, dimFrom+1};
          I = repmat((1:size(entity,1))',1,size(entity,2));
          obj.connectivity{dimFrom+1,dimTo+1} = spLeftShiftNonZero(sparse(entity, I, I))';
        end
        R = obj.connectivity{dimFrom+1,dimTo+1};
      end
    end
    function R = getNodePatch(obj, varargin) % [boundaryFlag]
      segm = obj.getEntity(1);
      if ~isempty(varargin) % on boundary
        iB = obj.isBoundary();
        segm = segm(iB,:);  
      end
      R = sparse(segm(:,[1 2]),segm(:,[2 1]),segm(:,[2 1]));
      R = spLeftShiftNonZero(R)';
    end
    function R = getProjector(obj) %#ok<*MANU>
      R = [];
    end
    function [R,cnt] = color(obj, dim)
      if dim==obj.dimP
        N = obj.getNumber(0);
        adjacency = obj.getEntity(obj.dimP);
        nColMax = size(adjacency,2)*max(accumarray(reshape(adjacency,[],1), 1, [N 1]));
      else
        N = obj.getNumber(obj.dimP);
        adjacency = obj.getConnect(dim, obj.dimP);
        nColMax = size(adjacency,2)*obj.nESub(dim+1);
      end
      % --> to MEX
%       color = zeros(N, nColMax);
%       [nA,nV] = size(adjacency);
%       R = zeros(nA,1);
%       for k = 1:nA
%         for idx = 1:nColMax
%           isFree = 0;
%           for l = 1:nV
%             if adjacency(k,l)==0, continue; end
%             isFree = isFree + color(adjacency(k,l), idx);
%           end
%           if isFree==0, break, end
%         end
%         R(k) = idx;
%         for l = 1:nV
%           if adjacency(k,l)==0, continue; end
%           color(adjacency(k,l), idx) = true;
%         end
%       end
      %
      R = color(adjacency, N, nColMax)+1;
      cnt = accumarray(R,1);
      [~,R] = sort(R);
    end
  end
  methods
    function adaptiveRefine(obj, I) %#ok<INUSD>
      error('adaptiveRefine() not yet implemented!');
    end
    function coarsen(obj, I) %#ok<INUSD>
      error('coarsen() not yet implemented!');
    end
  end
end