classdef GlobalSearcher < SOFE
  properties
    mesh
    nodes % nExnN
    dim
    bgMesh
    diam
    NVec
    nBlockGS = 1;
  end
  methods % constructor
    function obj = GlobalSearcher(mesh)
      obj.mesh = mesh;
      obj.dim = mesh.topology.dimP;
      obj.update();
    end
    function R = getBlock(obj, varargin) % [I]
      nE = obj.mesh.topology.getNumber('0');
      if obj.nBlockGS>nE, error('!Number of blocks exceeds number of elements!'); end
      R = unique(floor(linspace(0,nE,obj.nBlockGS+1)));
      R = [R(1:end-1)+1; R(2:end)];
      if nargin > 1
        R = R(:,varargin{1});
      end
    end
    function update(obj, varargin) % [nodes]
      if nargin < 2
        N = obj.mesh.nodes;
        elem = obj.mesh.topology.getEntity('0'); % nExnV
        I = (elem>0);
        obj.nodes = nan([numel(elem) obj.dim]);
        obj.nodes(I,:) = N(elem(I),:);
        obj.nodes = reshape(obj.nodes,size(elem,1),[],obj.dim); % nExnPxnW
      else
        obj.nodes = varargin{1}; % nExnPxnW
      end
      obj.diam(:,1) = min(reshape(obj.nodes,[],obj.dim)); % nWx2
      obj.diam(:,2) = max(reshape(obj.nodes,[],obj.dim)); % nWx2
      range = diff(obj.diam');
      obj.NVec = 5*obj.mesh.topology.getNumber('0')/2^obj.dim; % number of bins
      if obj.dim == 2
        obj.NVec = [obj.NVec*range(1)/range(2) ...
                    obj.NVec*range(2)/range(1)].^(1/2);
      elseif obj.dim == 3
        obj.NVec = [obj.NVec*range(1)^2/range(2)/range(3) ...
                    obj.NVec*range(2)^2/range(3)/range(1) ...
                    obj.NVec*range(3)^2/range(1)/range(2)].^(1/3);
      end
      obj.NVec = ceil(obj.NVec);
      obj.bgMesh = cell(obj.nBlockGS,1);
      repeat = true;
      while repeat
        repeat = false;
        for k = 1:obj.nBlockGS
          idx = obj.getBlock(k);
          obj.bgMesh{k} = obj.buildBackgroundMesh(idx(1):idx(2));
          if numel(obj.bgMesh{k}) == 1 && obj.bgMesh{k} < 0
            repeat = true;
            break
          end
          if k>1
            fprintf(repmat('\b',1,length(s)));
          end
          s = sprintf('Build background mesh %d / %d', k, obj.nBlockGS);
          fprintf(s);
        end
      end
      fprintf('\n');
    end
    function R = buildBackgroundMesh(obj, I)
      N = obj.nodes(I,:,:);
      bins = nan(size(N)); % nExnVxnW
      idx = ~isnan(N);
      bins(idx) = obj.getBin(reshape(N(idx),[],obj.dim));
      %
      minBins = permute(min(bins,[],2), [1 3 2]); % nExnW
      maxBins = permute(max(bins,[],2), [1 3 2]); % nExnW
      %
      nW = size(minBins, 2);
      XX = cell(nW, 1); iX = cell(nW, 1);
      for d = 1:nW
        XX{d} = reshapeTop(1-minBins(:,d)+maxBins(:,d)); % n{CoveredBins}_dxnE
        iX{d} = XX{d}==0;
        XX{d} = bsxfun(@plus, XX{d}, minBins(:,d)'-1);
        XX{d}(iX{d}) = 0;
      end
      %
      if nW == 1
        IDXE = [reshape(repmat(XX{1},1,1),[],1) ...
                reshape(kron((1:size(XX{1},2))', ones(size(XX{1},1),1)),[],1)];
      elseif nW == 2
        IDXE = [reshape(repmat(XX{1}, size(XX{2},1),1),[],1) ...
                reshape(kron(XX{2}, ones(size(XX{1},1),1)),[],1) ...
                reshape(kron((1:size(XX{1},2))', ones(size(XX{1},1)*size(XX{2},1),1)),[],1)];
      elseif nW == 3
        nComp = 4*numel(XX{1})*size(XX{2},1)*size(XX{3},1);
        if(nComp > 5e7)
          obj.nBlockGS = ceil(nComp/5e7);
          R = -1;
          return;
        end
        IDXE = zeros(numel(XX{1})*size(XX{2},1)*size(XX{3},1),4);
        IDXE(:,1) = reshape(repmat(XX{1}, size(XX{2},1)*size(XX{3},1),1),[],1);
        IDXE(:,2) = reshape(repmat(kron(XX{2}, ones(size(XX{1},1),1)), size(XX{3},1), 1),[],1);
        IDXE(:,3) = reshape(kron(XX{3}, ones(size(XX{1},1)*size(XX{2},1),1)), [], 1);
        IDXE(:,4) = reshape(kron((1:size(XX{1},2))', ones(size(XX{1},1)*size(XX{2},1)*size(XX{3},1),1)),[],1);
      end
      IDXE = IDXE(prod(IDXE(:,1:nW),2)>0,:); % (n{CoveredBins}_d*nE)x(nW+1)
      L = obj.getLinearIndex(IDXE(:,1:nW));
      IDXE = [L IDXE(:,nW+1)]; % (n{CoveredBins}_d*nE)x2
      [IDXE(:,1), pVec] = sort(IDXE(:,1));
      IDXE(:,2) = IDXE(pVec,2);
      if ~exist('OCTAVE_VERSION', 'builtin')
        [~, cnt] = unique(IDXE(:,1),'legacy');
      else
        [~, cnt] = unique(IDXE(:,1));
      end
      cnt = reshapeTop(diff([0;cnt]));
      cnt = cnt(cnt>0);
      R = sparse(IDXE(:,1), cnt, IDXE(:,2), prod(obj.NVec), max(cnt));
    end
  end
  methods % global search
    function R = findCandidates(obj, points)
      if obj.dim>1
        if isempty(obj.diam), obj.update(); end
        [~,L] = obj.getBin(points);
        I = L>0 & L<=prod(obj.NVec);
        R = cell(1,obj.nBlockGS);
        for k = 1:obj.nBlockGS
          idx = obj.getBlock(k);
          R{k} = zeros(numel(L), size(obj.bgMesh{k},2));
          R{k}(I,:) = obj.bgMesh{k}(L(I),:);
          R{k}(R{k}>0) = R{k}(R{k}>0) + (idx(1)-1);
        end
        R = cell2mat(R)';
        colSum = sum(R>0);
        R = reshapeTop(colSum, R(R>0))';
      else
        [pointsSorted, pid] = sort(points);
        R = zeros(length(points),1);
        elem = obj.mesh.topology.getEntity(1);
        [leftNodesSorted, nlid] = sort(obj.mesh.nodes(elem(:,1)));
        if (leftNodesSorted(1)<pointsSorted(1))
          j = 1;
        else
          j = 2;
        end
        for i=1:length(points)
          while (j<=length(leftNodesSorted)) && (leftNodesSorted(j)<pointsSorted(i))
            j=j+1;
          end
          R(pid(i)) = nlid(j-1);
        end
      end
    end
    function [R, L] = getBin(obj, points)
      sPoints = bsxfun(@rdivide, bsxfun(@minus, points, obj.diam(:,1)'), diff(obj.diam')); % nPx3
      sPoints = max(0,min(1-eps,sPoints));
      R = floor(bsxfun(@times, obj.NVec, sPoints)) + 1;
      L = obj.getLinearIndex(R);
    end
    function R = getLinearIndex(obj, idx)
      nW = obj.dim;
      correct = [0, cumprod(obj.NVec)];
      tmp = [1 cumprod(obj.NVec)];
      R = sum(bsxfun(@times, idx, tmp(1:nW)), 2) - sum(correct(1:nW));
    end
  end
end