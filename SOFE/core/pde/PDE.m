classdef PDE < SOFE
  properties
    nEq, nOp
    list, lhs, rhs
    %
    stiffMat, loadVec, shift
    freeI, freeJ
    createSys = true;
    %
    mesh
    fesTest, fesTrial
    I,J, nDoF
    %
    time, state
    nArgIn
    stateChanged
  end
  methods % constructor & more
    function obj = PDE(list, lhs, rhs)
      obj.list = list; obj.nOp = numel(list);
      obj.lhs = lhs; obj.rhs = rhs; obj.nEq = numel(rhs.sys);
      obj.fesTest = cell(obj.nEq, 1); obj.fesTrial = cell(obj.nEq, 1);
      for k = 1:numel(list)
        obj.list{k}.register(obj);
      end
      for i = 1:obj.nEq
        for j = 1:obj.nEq
          if ~isempty(obj.lhs.sys{i,j})
            fes{2} = obj.list{obj.lhs.sys{i,j}{1}}.fesTest;
            fes{1} = obj.list{obj.lhs.sys{i,j}{1}}.fesTrial;
            try adj = obj.lhs.adj{i,j}{1}; catch, adj = false; end
            if adj; fes = fes([2 1]); end
            if isempty(obj.fesTrial{j})
              obj.fesTrial{j} = fes{1};
            end
            if isempty(obj.fesTest{i})
              obj.fesTest{i} = fes{2};
            end
          end
        end
      end
      obj.nArgIn.coeff =  max(cellfun(@(op)nargin(op.dataCache),obj.list));
      obj.nArgIn.shift =  max(cellfun(@(fes)fes.narginShift,obj.fesTrial));
      obj.mesh = obj.fesTrial{1}.mesh;
      obj.notify('init');
    end
    function notify(obj, varargin) % [message]
      [nTest, nTrial] = obj.getNDoF();
      nTest = cumsum(nTest);
      obj.I = [[1;nTest(1:end-1)+1], nTest];
      nTrial = cumsum(nTrial);
      obj.J = [[1;nTrial(1:end-1)+1], nTrial];
      %
      obj.nDoF = nTrial(obj.nEq);
      obj.stiffMat = [];
      obj.loadVec = []; obj.shift = [];
      obj.freeI = []; obj.freeJ = [];
      obj.setState(0.0, zeros(obj.nDoF,1));
    end
    function setState(obj, t, varargin) % [state]
      obj.time = t;
      obj.stateChanged = true;
      if ~isempty(varargin)
        obj.state = varargin{1};
      end
      if obj.nArgIn.shift>1, obj.shift = []; end
    end
    function setMatrixFree(obj, hasCoeff)
      obj.createSys = 0;
      for k = 1:numel(obj.list)
        obj.list{k}.matrixFree = 1;
        if ~isa(obj.list{k}, 'Operator')
          obj.list{k}.hasCoeff = 0;
        else
          obj.list{k}.hasCoeff = hasCoeff;
        end
      end
      obj.assemble();
    end
  end
  methods
    function assemble(obj)
      if obj.stateChanged && (isempty(obj.loadVec) || isempty(obj.stiffMat) || (obj.nArgIn.coeff > 1))
        obj.stateChanged = false;
        obj.stiffMat = []; obj.loadVec = [];
        for k = 1:obj.nOp
          obj.list{k}.notify('stateChanged', obj.time);
          obj.list{k}.assemble();
        end
        obj.createSystem();
      end
    end
    function createSystem(obj)
      if nnz(obj.stiffMat)>0
        fprintf('System already created\n');
        keyboard
        return;
      end
      if obj.createSys, obj.stiffMat = sparse(obj.nDoF, obj.nDoF); end
      obj.loadVec = zeros(obj.nDoF, 1);
      for i = 1:obj.nEq
        % rhs
        if ~isempty(obj.rhs.sys{i}) 
          for k = 1:numel(obj.rhs.sys{i})
            b = obj.list{obj.rhs.sys{i}{k}}.matrix;
            try b = obj.rhs.coeff{i}{k}*b; catch, end
            idx = obj.I(i,1):obj.I(i,2);
            obj.loadVec(idx,:) = obj.loadVec(idx,:) + b;
          end
        end
        % lhs
        if obj.createSys
          for j = 1:obj.nEq
            if ~isempty(obj.lhs.sys{i,j})
              for k = 1:numel(obj.lhs.sys{i,j})
                try cc = obj.lhs.coeff{i,j}{k}; catch, cc = 1.0; end
                try adj = obj.lhs.adj{i,j}{k}; catch, adj = 0; end
                if adj
                  obj.list{obj.lhs.sys{i,j}{k}}.matrix = obj.list{obj.lhs.sys{i,j}{k}}.matrix.';
                end
                obj.stiffMat = obj.stiffMat + ...
                       [sparse(obj.nDoF, obj.J(j,1)-1), ...
                       [sparse(obj.I(i,1)-1, obj.J(j,2)-obj.J(j,1)+1); ...
                        cc*obj.list{obj.lhs.sys{i,j}{k}}.matrix; ...
                        sparse(obj.nDoF-obj.I(i,2), obj.J(j,2)-obj.J(j,1)+1)], ...
                        sparse(obj.nDoF, obj.nDoF-obj.J(j,2))];
                if adj
                  obj.list{obj.lhs.sys{i,j}{k}}.matrix = obj.list{obj.lhs.sys{i,j}{k}}.matrix.';
                end
              end
            end
          end
        end
      end
    end
    function R = applySystem(obj, x, varargin) % [isFree]
      if ~isempty(varargin)
        [fI, fJ] = obj.getFreeDoFs();
      else
        fI = ':'; fJ = ':';
      end
      %
      R = zeros(obj.nDoF, 1);
      X = zeros(obj.nDoF, 1);
      X(fJ) = x;
      if ~isempty(obj.stiffMat)
        R = obj.stiffMat*X;
      else
        for i = 1:obj.nEq
          for j = 1:obj.nEq
            if ~isempty(obj.lhs.sys{i,j})
              for k = 1:numel(obj.lhs.sys{i,j})
                try adj = obj.lhs.adj{i,j}{k}; catch, adj = 0; end
                if adj
                  dR = (X(obj.J(j,1):obj.J(j,2))'*obj.list{obj.lhs.sys{i,j}{k}}.matrix)';
                else
%                   dR = obj.list{obj.lhs.sys{i,j}{k}}.matrix*X(obj.J(j,1):obj.J(j,2));
                  dR = obj.list{obj.lhs.sys{i,j}{k}}.apply(X(obj.J(j,1):obj.J(j,2)));
                end
                try dR = obj.lhs.coeff{i,j}{k}*dR; catch, end
                R(obj.I(i,1):obj.I(i,2)) = R(obj.I(i,1):obj.I(i,2)) + dR;
              end
            end
          end
        end
      end
      R = R(fI);
    end
    function R = applyBlock(obj, kl, x, varargin) % [isFree, transposeFlag]
      k = kl(1); l = kl(2);
      if ~isempty(varargin) && ~isempty(varargin{1})
        [fI, ~] = obj.getFreeDoFs(k);
        [~, fJ] = obj.getFreeDoFs(l);
        if numel(varargin)>1
          tmp = fI; fI = fJ; fJ = tmp;
        end
      else
        fI = ':'; fJ = ':';
      end
      %
      R = zeros(size(fI));
      X = zeros(size(fJ));
      X(fJ) = x;
      for m = 1:numel(obj.lhs.sys{k,l})
        try adj = obj.lhs.adj{k,l}{m}; catch, adj = 0; end
        if numel(varargin)>1, adj = ~adj; end
        if adj
          dR = (X'*obj.list{obj.lhs.sys{k,l}{m}}.matrix)';
        else
          dR = obj.list{obj.lhs.sys{k,l}{m}}.matrix*X;
        end
        try dR = obj.lhs.coeff{k,l}{m}*dR; catch, end
        R = R + dR;
      end
      R = R(fI);
    end
    function R = evalState(obj, k) % [k]
      R.U = cell(obj.nEq, 1); % {nEq}xnExnPxnC
      R.dU = cell(obj.nEq, 1); % {nEq}xnExnPxnCxnD
      if obj.nArgIn.coeff > 2
        for j = 1:obj.nEq
          U = obj.getStateCmp(j);
          R.U{j} = obj.fesTrial{j}.evalDoFVector(U,[],0,0, {k});
          if obj.nArgIn.coeff > 3
            R.dU{j} = obj.fesTrial{j}.evalDoFVector(U,[],0,1, {k});
          end
        end
      end
    end
    function R = getShift(obj, varargin) % [idx]
      if isempty(obj.shift)
        obj.shift = cellfun(@(fes)fes.getShift(obj.time),obj.fesTrial,'UniformOutput',0);
      end
      if ~isempty(varargin)
        R = obj.shift{varargin{1}};
      else
        R = cell2mat(obj.shift);
      end
    end
    function R = getStateCmp(obj, varargin) % [idx]
      if isempty(varargin)
        R = obj.state;
      else
        R = obj.state(obj.J(varargin{1},1):obj.J(varargin{1},2));
      end
    end
    function [RI, RJ] = getFreeDoFs(obj, varargin) % [idx]
      if isempty(obj.freeI)
        obj.freeI = cell2mat(cellfun(@(fes)fes.getFreeDoFs(),obj.fesTest,'UniformOutput',0));
        obj.freeJ = cell2mat(cellfun(@(fes)fes.getFreeDoFs(),obj.fesTrial,'UniformOutput',0));
      end
      if ~isempty(varargin)
        idx = varargin{1};
        RI = obj.freeI(obj.I(idx,1):obj.I(idx,2));
        RJ = obj.freeJ(obj.J(idx,1):obj.J(idx,2));
      else
        RI = obj.freeI;
        RJ = obj.freeJ;
      end
    end
    function [RI, RJ] = getNDoF(obj, varargin) % [idx]
      RI = cellfun(@(fes)fes.getNDoF,obj.fesTest);
      RJ = cellfun(@(fes)fes.getNDoF,obj.fesTrial);
      if ~isempty(varargin)
        RI = RI(varargin{1});
        RJ = RJ(varargin{1});
      end
    end
  end
end