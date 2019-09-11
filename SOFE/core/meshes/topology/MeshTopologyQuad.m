classdef MeshTopologyQuad < MeshTopology
  methods % constructor
    function obj = MeshTopologyQuad(elem)
      obj = obj@MeshTopology(2);
      obj.update(elem);
      obj.isSimplex = 0;
      obj.nESub = [4 4 1];
      obj.nO = 2;
    end
    function update(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{2,1} = [elem(:,[1,2]); elem(:,[3,4]); elem(:,[1,3]); elem(:,[2,4])];
      [obj.connectivity{2,1}, ~, e2F] = unique(sort(obj.connectivity{2,1},2),'rows'); 
      obj.connectivity{3,2} = reshape(e2F, size(elem,1), []);
      %
      obj.connectivity{1,1} = (1:max(obj.connectivity{2,1}(:)))';
      obj.connectivity{2,2} = (1:size(obj.connectivity{2,1},1))';
      obj.connectivity{3,3} = (1:size(obj.connectivity{3,1},1))';
    end
  end
  methods % connectivity information   
    function R = getOrientation(obj, dim, d, varargin) % [I]
      R = [];
      if dim==2 && d == 1
        e = obj.getEntity('0', varargin{:});
        R = ones(size(e));
        R(e(:,1)>e(:,2),1) = 2;
        R(e(:,3)>e(:,4),2) = 2;
        R(e(:,1)>e(:,3),3) = 2;
        R(e(:,2)>e(:,4),4) = 2;
      end
    end
    function R = getNormalOrientation(obj, varargin) % [I]
      R = 2*mod(obj.getOrientation(2, 1, varargin{:}),2)-1; % nExnF
      R(:,[2 3]) = -R(:,[2 3]);
    end
  end
  methods % refinement
    function R = uniformRefine_(obj)
      fc = obj.getEntity(1);
      el = obj.getEntity(2);
      nN = obj.getNumber(0); nF = obj.getNumber(1); nE = obj.getNumber(2);
      R = [speye(nN); sparse(repmat((1:nF)',1,2), fc, 0.5); ...
                    sparse(repmat((1:nE)',1,4), el, 0.25)];
      newIndicesF = nN + (1:nF);
      newIndicesE = nN + nF + (1:nE)';
      el = [el newIndicesF(obj.connectivity{3,2}) newIndicesE];
      el = [el(:,[1 5 7 9]);el(:,[5 2 9 8]);el(:,[7 9 3 6]);el(:,[9 8 6 4])];
      obj.update(el);
    end
    function P = uniformRefine(obj)
      fc = obj.getEntity(1); el = obj.getEntity(2);
      e2F = obj.connectivity{3,2}; oo = obj.getOrientation(2,1)==2;
      nN = obj.getNumber(0); nF = obj.getNumber(1); nE = obj.getNumber(2);
      fRange = (1:nF)'; eRange = (1:nE)';
      %
      P = [speye(nN); fsparse(repmat((1:nF)',1,2), fc, 0.5); fsparse(repmat((1:nE)',1,4), el, 0.25)];
      %
      newIndicesF = nN + fRange; newIndicesE = nN + nF + eRange;
      el = [el reshape(newIndicesF(obj.connectivity{3,2}),size(obj.connectivity{3,2})) newIndicesE];
      el = [el(:,[1 5 7 9]);el(:,[5 2 9 8]);el(:,[7 9 3 6]);el(:,[9 8 6 4])];
      fc = [fc(:,1) newIndicesF; fc(:,2) newIndicesF; ...
            nN+e2F(:,1) newIndicesE; nN+e2F(:,2) newIndicesE; nN+e2F(:,3) newIndicesE; nN+e2F(:,4) newIndicesE];
      e2F = [[nF*oo(:,1)+e2F(:,1), 2*(nF+nE)+eRange, nF*oo(:,3)+e2F(:,3), 2*nF+eRange]; ...
             [nF*~oo(:,1)+e2F(:,1), (2*nF+3*nE)+eRange, 2*nF+eRange, nF*oo(:,4)+e2F(:,4)]; ...
             [2*(nF+nE)+eRange, nF*oo(:,2)+e2F(:,2), nF*~oo(:,3)+e2F(:,3), (2*nF+nE)+eRange]; ...
             [(2*nF+3*nE)+eRange, nF*~oo(:,2)+e2F(:,2), (2*nF+nE)+eRange, nF*~oo(:,4)+e2F(:,4)]];
      %
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{3,1} = el;
      obj.connectivity{2,1} = fc;
      obj.connectivity{3,2} = e2F;
      obj.connectivity{1,1} = (1:size(P,1))';
      obj.connectivity{2,2} = (1:size(fc,1))';
      obj.connectivity{3,3} = (1:size(el,1))';
    end
  end
  methods(Static = true)
    function R = isFeasible(points, varargin) % [tol
      if ~isempty(varargin), tol = varargin{1}; else, tol = 1e-12; end
      R = (points(:,1)>-tol & points(:,1)<1+tol & points(:,2)>-tol & points(:,2)<1+tol);
    end
    function R = getCenterLoc()
      R = [1 1]/2;
    end
    function R = upliftPoints(points, fLoc, orient)
      % complies standard orientation
      zz = zeros(size(points)); oo = ones(size(points));
      if orient==2
        points = 1-points;
      end
      switch fLoc
        case 1
          R = [points, zz];
        case 2
          R = [points, oo];
        case 3
          R = [zz, points];
        case 4
          R = [oo, points];
      end
    end
  end
end