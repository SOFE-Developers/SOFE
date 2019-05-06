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
    function [R, nType] = getNodePatch(obj, dim)
      if dim>0
        if isempty(obj.connectivity{1,dim+1}) || nargout>1
          entity = obj.getEntity(dim);
          idx = repmat((1:size(entity,1))', 1, size(entity,2));
          entity = entity(:);
          [~,I] = sort(entity);
          count = accumarray(entity,1);
          maxCount = max(count);
          upperTri = triu(repmat((1:maxCount)',1,maxCount));
          count = upperTri(:,count); count = count(:);
          count(count==0) = [];
          obj.connectivity{1,dim+1} = accumarray([entity(I), count], idx(I));
          if nargout>1
            idx = kron(ones(size(idx,1),1), (1:size(idx,2)));
            nType = accumarray([entity(I), count], idx(I));
          end
        end
        R = obj.connectivity{1,dim+1};
      else
        segm = obj.getEntity(1);
        if dim < 0 % on boundary
          iB = obj.isBoundary();
          segm = segm(iB,:);  
        end
        R = sparse(segm(:,[1 2]),segm(:,[2 1]),segm(:,[2 1]));
        R = spLeftShiftNonZero(R)';
%         if dim < 0 % on boundary
%           R = full(II(:,1:2));
%         else
%           R = full(II(:,1:max(sum(II>0,2))));
%         end
%         R(R(:,1)==0,:) = []; % deprecated
      end
    end
    function R = getNodePatch2(obj, dim)
      E = obj.getEntity(dim);
      I = repmat((1:size(E,1))',1,size(E,2));
      R = sparse(E, I, I);
      R = spLeftShiftNonZero(R)';
    end
    function R = getProjector(obj) %#ok<*MANU>
      R = [];
    end
    function [R,cnt] = getColoredElements(obj)
      nN = obj.getNumber(0);
      elem = obj.getEntity('0');
      degN = accumarray(reshape(elem,[],1), 1, [nN 1]);
      nColMax = size(elem,2)*max(degN);
      % --> to MEX
%       nodeColor = zeros(nN, nColMax);
%       R = zeros(size(elem,1),1);
%       for k = 1:size(elem,1)
%         for idx = 1:nColMax
%           isFree = 0;
%           for l = 1:size(elem,2)
%             isFree = isFree + nodeColor(elem(k,l), idx);
%           end
%           if isFree==0, break, end
%         end
%         R(k) = idx;
%         nodeColor(elem(k,:), idx) = true;
%       end
      %
      R = colorElements(elem, nN, nColMax)+1;
      cnt = accumarray(R,1);
      [~,R] = sort(R);
    end
  end
  methods
    function uniformRefineFast(obj, I) %#ok<*INUSD>
      error('uniformRefineFast() not yet implemented!');
    end
    function adaptiveRefine(obj, I)
      error('adaptiveRefine() not yet implemented!');
    end
    function coarsen(obj, I)
      error('coarsen() not yet implemented!');
    end
  end
end