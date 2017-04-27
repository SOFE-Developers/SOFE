classdef GlobalSearcher < SOFE
  properties
    topology
    bgMesh
    diam
    NVec
    nBlockGS = 1;
  end
  methods % constructor
    function obj = GlobalSearcher(topology)
      obj.topology = topology;
      fprintf('Set up GlobalSearcher ... \n');
      obj.notify();
      fprintf('DONE\n');
    end
    function R = getBlock(obj, varargin) % [I]
      nE = obj.topology.getNumber(obj.topology.dimP);
      if obj.nBlockGS>nE, error('!Number of blocks exceeds number of elements!'); end
      R = unique(floor(linspace(0,nE,obj.nBlockGS+1)));
      R = [R(1:end-1)+1; R(2:end)];
      if nargin > 1
        R = R(:,varargin{1});
      end
    end
    function notify(obj)
      dim = obj.topology.dimW;
      obj.diam(:,1) = min(obj.topology.nodes); % nWx2
      obj.diam(:,2) = max(obj.topology.nodes); % nWx2
      range = diff(obj.diam');
      obj.diam(:,1) = obj.diam(:,1) - 0.01*range'; % nWx2
      obj.diam(:,2) = obj.diam(:,2) + 0.01*range'; % nWx2
      obj.NVec = obj.topology.getNumber(obj.topology.dimP)/2^dim; % number of bins
      if dim == 2
        obj.NVec = [obj.NVec*range(1)/range(2) ...
                    obj.NVec*range(2)/range(1)].^(1/2);
      elseif dim == 3
        obj.NVec = [obj.NVec*range(1)^2/range(2)/range(3) ...
                    obj.NVec*range(2)^2/range(3)/range(1) ...
                    obj.NVec*range(3)^2/range(1)/range(2)].^(1/3);
      end
      obj.NVec = ceil(obj.NVec);
      obj.bgMesh = cell(obj.nBlockGS,1);
      for k = 1:obj.nBlockGS
        idx = obj.getBlock(k);
        obj.bgMesh{k} = obj.buildBackgroundMesh(idx(1):idx(2));
      end
    end
    function R = buildBackgroundMesh(obj, I)
      elem = obj.topology.getEntity(obj.topology.dimP, I); % nExnV
      bins = zeros([size(elem) size(obj.topology.nodes,2)]); % nExnVxnW
      nV = size(elem,2);
      for k = 1:nV
       bins(:,k,:) = obj.getBin(obj.topology.nodes(elem(:,k),:)); % nExnVxnW
      end
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
        IDXE = zeros(prod(size(XX{1}))*size(XX{2},1)*size(XX{3},1),4);
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
    end
    function [R, L] = getBin(obj, points)
      sPoints = bsxfun(@rdivide, bsxfun(@minus, points, obj.diam(:,1)'), diff(obj.diam')); % nPx3
      R = floor(bsxfun(@times, obj.NVec, sPoints)) + 1;
      L = obj.getLinearIndex(R);
    end
    function R = getLinearIndex(obj, idx)
      nW = obj.topology.dimW;
      correct = [0, cumprod(obj.NVec)];
      tmp = [1 cumprod(obj.NVec)];
      R = sum(bsxfun(@times, idx, tmp(1:nW)), 2) - sum(correct(1:nW));
    end
  end
end