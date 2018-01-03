classdef PDE < SOFE
  properties
    nEq, nOp
    list, lhs, rhs
    %
    stiffMat, loadVec
    createSys = true;
    %
    mesh
    fesTest, fesTrial
    I,J, nDoF
    %
    time, state
    narginData
    stateChanged
  end
  methods % constructor & more
    function obj = PDE(list, lhs, rhs)
      obj.list = list; obj.nOp = numel(list);
      obj.lhs = lhs; obj.rhs = rhs; obj.nEq = numel(rhs.sys);
      obj.fesTest = cell(obj.nEq, 1); obj.fesTrial = cell(obj.nEq, 1);
      obj.narginData =  max(cellfun(@(op)nargin(op.dataCache),obj.list)); % maximal nargin of data(x,t,u,d)
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
      obj.notify();
      obj.mesh = obj.fesTrial{1}.mesh;
      obj.setState(0.0, zeros(obj.nDoF,1));
    end
    function notify(obj)
      [nTest, nTrial] = obj.getNDoF();
      nTest = cumsum(nTest);
      obj.I = [[1;nTest(1:end-1)+1], nTest];
      nTrial = cumsum(nTrial);
      obj.J = [[1;nTrial(1:end-1)+1], nTrial];
      %
      obj.nDoF = nTrial(obj.nEq);
      obj.stiffMat = []; obj.loadVec = [];
    end
    function setState(obj, t, varargin) % [state]
      obj.time = t;
      obj.stateChanged = true;
      if ~isempty(varargin)
        obj.state = varargin{1};
      end
    end
  end
  methods
    function assemble(obj)
      if obj.stateChanged && (isempty(obj.stiffMat) || (obj.narginData > 1))
        obj.stateChanged = false;
        obj.stiffMat = []; obj.loadVec = [];
        for k = 1:obj.nOp
          obj.list{k}.notify(obj.time);
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
      obj.stiffMat = sparse(obj.nDoF, obj.nDoF);
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
    function R = applySystem(obj, x, varargin) % [onFreeDoFs or {freeI, freeJ}]
      if isempty(varargin)
        freeI = ':'; freeJ = ':';
      else
        freeI = varargin{1}; freeJ = varargin{2};
      end
      R = zeros(obj.nDoF, 1);
      X = zeros(obj.nDoF, 1);
      X(freeJ) = x;
      for i = 1:obj.nEq
        for j = 1:obj.nEq
          if ~isempty(obj.lhs.sys{i,j})
            for k = 1:numel(obj.lhs.sys{i,j})
              try adj = obj.lhs.adj{i,j}{k}; catch, adj = 0; end
              if adj
                dR = (X(obj.J(j,1):obj.J(j,2))'*obj.list{obj.lhs.sys{i,j}{k}}.matrix)';
              else
                %dR = obj.list{obj.lhs.sys{i,j}{k}}.matrix*X(obj.J(j,1):obj.J(j,2));
                dR = obj.list{obj.lhs.sys{i,j}{k}}.apply(X(obj.J(j,1):obj.J(j,2)));
              end
              try dR = obj.lhs.coeff{i,j}{k}*dR; catch, end
              R(obj.I(i,1):obj.I(i,2)) = R(obj.I(i,1):obj.I(i,2)) + dR;
            end
          end
        end
      end
      R = R(freeI);
    end
    function R = evalState(obj, k) % [k]
      R.U = cell(obj.nEq, 1); % {nEq}xnExnPxnC
      R.dU = cell(obj.nEq, 1); % {nEq}xnExnPxnCxnD
      if obj.narginData > 2
        for j = 1:obj.nEq
          U = obj.state(obj.J(j,1):obj.J(j,2));
          R.U{j} = obj.fesTrial{j}.evalDoFVector(U,[],0,0, {k});
          if obj.narginData > 3
            R.dU{j} = obj.fesTrial{j}.evalDoFVector(U,[],0,1, {k});
          end
        end
      end
    end
    function R = getShift(obj, varargin) % [idx]
      R = cellfun(@(fes)fes.getShift(obj.time),obj.fesTrial,'UniformOutput',0);
      if ~isempty(varargin)
        R = R{varargin{1}};
      else
        R = cell2mat(R);
      end
    end
    function R = getState(obj, varargin) % [idx]
      if isempty(varargin)
        R = obj.state;
      else
        R = obj.state(obj.J(varargin{1},1):obj.J(varargin{1},2));
      end
    end
    function [RI, RJ] = getFreeDoFs(obj, varargin) % [idx]
      RI = cellfun(@(fes)fes.getFreeDoFs(),obj.fesTest,'UniformOutput',0);
      RJ = cellfun(@(fes)fes.getFreeDoFs(),obj.fesTrial,'UniformOutput',0);
      if ~isempty(varargin)
        RI = RI{varargin{1}};
        RJ = RJ{varargin{1}};
      else
        RI = cell2mat(RI);
        RJ = cell2mat(RJ);
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