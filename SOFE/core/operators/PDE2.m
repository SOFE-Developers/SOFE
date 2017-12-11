classdef PDE2 < SOFE
  properties
    nEq, nOp
    list, lhs, rhs
    %
    stiffMat, loadVec, solution
    shift, fDoFsTrial, fDoFsTest
    %
    mesh
    fesTest,fesTrial
    I,J, nDoF
    %
    time, state, dState
    narginData
    stateChanged
    %
    solver = DirectSolver([]);
  end
  methods % constructor & more
    function obj = PDE2(list, lhs, rhs)
      obj.list = list;
      obj.lhs = lhs;
      obj.rhs = rhs;
      obj.nEq = numel(rhs.sys);
      obj.nOp = numel(list);
      obj.fesTrial = cell(obj.nEq, 1);
      obj.fesTest = cell(obj.nEq, 1);
      obj.narginData = 1; % maximal nargin of data(x,t,u,d)
      for k = 1:obj.nOp
        obj.list{k}.pde = obj;
        obj.narginData = max(obj.narginData, nargin(obj.list{k}.dataCache));
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
              obj.fesTrial{j}.register(obj);
            end
            if isempty(obj.fesTest{i})
              obj.fesTest{i} = fes{2};
              obj.fesTest{i}.register(obj);
            end
          end
        end
      end
      obj.mesh = obj.fesTrial{1}.mesh;
      obj.setIndices();
      obj.setTime(0.0);
      obj.setState();
    end
    function setIndices(obj)
      dimTest = zeros(obj.nEq, 1); dimTrial = zeros(obj.nEq, 1);
      for i = 1:obj.nEq
        dimTest(i) = obj.fesTest{i}.getNDoF();
        dimTrial(i) = obj.fesTrial{i}.getNDoF();
      end
      offsetTrial = [0; cumsum(dimTrial)]; offsetTest = [0; cumsum(dimTest)];
      obj.J = cell(obj.nEq, 1); obj.I = cell(obj.nEq, 1);
      for i = 1:obj.nEq
        obj.J{i} = offsetTrial(i) + [1 dimTrial(i)];
        obj.I{i} = offsetTest(i) + [1 dimTest(i)];
      end
      obj.nDoF = obj.I{obj.nEq}(2);
    end
    function notify(obj)
      obj.stiffMat = []; obj.loadVec = []; obj.solution = [];
      obj.shift = []; obj.fDoFsTest = []; obj.fDoFsTrial = [];
      obj.setIndices(); obj.setState();
    end
    function setTime(obj, newTime)
      obj.stateChanged = true;
      obj.time = newTime;
      if obj.narginData>1, obj.stiffMat = []; obj.loadVec = []; end
      for k = 1:obj.nEq
        if ~isempty(obj.fesTrial{k}.shift) && nargin(obj.fesTrial{k}.shift) > 1
          obj.shift = [];
        end
        if ~isempty(obj.fesTrial{k}.fixB) && nargin(obj.fesTrial{k}.fixB) > 1
          obj.fDoFsTrial = [];
        end
        if ~isempty(obj.fesTest{k}.fixB) && nargin(obj.fesTest{k}.fixB) > 1
          obj.fDoFsTest = [];
        end
      end
    end
    function setState(obj, varargin) % [U]
      obj.stateChanged = true;
      if obj.narginData < 3, return; end
      obj.stiffMat = []; obj.loadVec = [];
      obj.state = cell(obj.nEq, 1); % {nEq}xnExnP
      for j = 1:obj.nEq
        [~, w] = obj.fesTrial{j}.getQuadData(0);
        nD = obj.fesTrial{j}.element.dimension;
        nC = size(obj.fesTrial{j}.element.evalBasis(zeros(1,nD),0),3);
        nE = obj.mesh.topology.getNumber('0');
        obj.state{j} = zeros(nE, numel(w), nC);
        if nargin < 2, continue; end
        U = varargin{1}(obj.J{j}(1):obj.J{j}(2));
        obj.state{j} = obj.fesTrial{j}.evalDoFVector(U, [], 0, 0);
      end
      if obj.narginData < 4, return; end
      obj.dState = cell(obj.nEq, 1); % {nEq}xnExnP
      for j = 1:obj.nEq
        [~, w] = obj.fesTrial{j}.getQuadData(0);
        nD = obj.fesTrial{j}.element.dimension;
        nC = size(obj.fesTrial{j}.element.evalBasis(zeros(1,nD),0),3);
        nE = obj.mesh.topology.getNumber('0');
        obj.dState{j} = zeros(nE, numel(w), nC, nD);
        if nargin < 2, continue; end
        U = varargin{1}(obj.J{j}(1):obj.J{j}(2));
        obj.dState{j} = obj.fesTrial{j}.evalDoFVector(U, [], 0, 1);
      end
    end
  end
  methods % assemble & solve
    function compute(obj)
      t = tic; obj.output('Begin assemble ...', 1);
      obj.assemble();
      fprintf('%d DoFs\n', sum(obj.fDoFsTrial));
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);      
      t = tic; obj.output('Begin solve ...', 1);
      obj.solve();
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
    function assemble(obj)
      obj.stateChanged = false;
      if ~isempty(obj.stiffMat) && (~obj.stateChanged || obj.narginData < 2)
        return
      end
      % operators
      for k = 1:obj.nOp
        obj.list{k}.notify(obj.time);
        obj.list{k}.assemble();
      end
      % freeDofs
      if isempty(obj.shift)
        obj.fDoFsTest = true(obj.I{obj.nEq}(2), 1);
        obj.fDoFsTrial = true(obj.J{obj.nEq}(2), 1);
        obj.shift = zeros(obj.J{obj.nEq}(2), 1);      
        for j = 1:obj.nEq
          obj.shift(obj.J{j}(1):obj.J{j}(2),1) = obj.fesTrial{j}.getShiftVector(obj.time);
          obj.fDoFsTest(obj.I{j}(1):obj.I{j}(2),1) = obj.fesTest{j}.getFreeDoFs();
          obj.fDoFsTrial(obj.J{j}(1):obj.J{j}(2),1) = obj.fesTrial{j}.getFreeDoFs();
        end
      end
      obj.createSystem();
    end
    function createSystem(obj)
      if ~isempty(obj.stiffMat)
        fprintf('System already created\n');
        keyboard
        return;
      end
      obj.stiffMat = sparse(obj.I{obj.nEq}(2), obj.J{obj.nEq}(2));
      obj.loadVec = zeros(obj.I{obj.nEq}(2), 1);
      for i = 1:obj.nEq
        % lhs
        for j = 1:obj.nEq
          if ~isempty(obj.lhs.sys{i,j})
            for k = 1:numel(obj.lhs.sys{i,j})
              try cc = obj.lhs.coeff{i,j}{k}; catch, cc = 1; end
              try adj = obj.lhs.adj{i,j}{k}; catch, adj = 0; end
              if adj
                obj.list{obj.lhs.sys{i,j}{k}}.matrix = obj.list{obj.lhs.sys{i,j}{k}}.matrix.';
              end
              obj.stiffMat = obj.stiffMat + ...
                     [sparse(obj.I{obj.nEq}(2), obj.J{j}(1)-1), ...
                     [sparse(obj.I{i}(1)-1, obj.J{j}(2)-obj.J{j}(1)+1); ...
                      cc*obj.list{obj.lhs.sys{i,j}{k}}.matrix; ...
                      sparse(obj.I{obj.nEq}(2)-obj.I{i}(2), obj.J{j}(2)-obj.J{j}(1)+1)], ...
                      sparse(obj.I{obj.nEq}(2), obj.J{obj.nEq}(2)-obj.J{j}(2))];
              if adj
                obj.list{obj.lhs.sys{i,j}{k}}.matrix = obj.list{obj.lhs.sys{i,j}{k}}.matrix.';
              end
            end
          end
        end
        % rhs
        if ~isempty(obj.rhs.sys{i})
          for k = 1:numel(obj.rhs.sys{i})
            b = obj.list{obj.rhs.sys{i}{k}}.vector;
            try cc = obj.lhs.coeff{i,j}{k}; catch, cc = 1; end
            idx = obj.I{i}(1):obj.I{i}(2);
            obj.loadVec(idx,:) = obj.loadVec(idx,:) + cc*b;
          end
        end
      end
    end
    function R = applySystem(obj, x, varargin) % [freeI, freeJ]
      if isempty(varargin)
        freeI = obj.fDoFsTest; freeJ = obj.fDoFsTrial;
      else
        freeI = varargin{1}; freeJ = varargin{2};
      end
      R = zeros(size(x));
      for i = 1:obj.nEq
        for j = 1:obj.nEq
          if ~isempty(obj.lhs.sys{i,j})
            for k = 1:numel(obj.lhs.sys{i,j})
              II = obj.I{i}(1):obj.I{i}(2); JJ = obj.J{j}(1):obj.J{j}(2);
              ii = freeI(II)>0; jj = freeJ(JJ)>0;
              II(~freeI(II)) = []; JJ(~freeJ(JJ)) = [];
              R(II) = R(II) + obj.lhs{i,j}{k}.matrix(ii,jj)*x(JJ);
            end
          end
        end
      end
    end
    %
    function solve(obj)
      obj.solution = zeros(size(obj.loadVec));
      for k = 1:size(obj.loadVec,2)
        b = obj.loadVec(:,k);
        b = b - obj.stiffMat*obj.shift; 
        obj.solution(~obj.fDoFsTrial, k) = obj.shift(~obj.fDoFsTrial);
        A = obj.stiffMat(obj.fDoFsTest, obj.fDoFsTrial);
        x = obj.solver.solve(A, b(obj.fDoFsTest));
        obj.solution(obj.fDoFsTrial, k) = obj.shift(obj.fDoFsTrial) + x;
      end
    end
  end
  methods % access
    function R = getTrialSpace(obj, j)
      R = obj.lhs{1,j}{1}.fesTrial;
    end
    function R = getTestSpace(obj, i)
      R = obj.lhs{i,1}{1}.fesTest;
    end
    function R = getSolution(obj, idx)
      R = obj.solution(obj.J{idx}(1):obj.J{idx}(2));
    end
    function setSolution(obj, sol, idx)
      obj.solution(obj.J{idx}(1):obj.J{idx}(2)) = sol;
    end
    function R = getStiffnessBlock(obj, idxI, idxJ)
      R = obj.stiffMat(obj.I{idxI},obj.J{idxJ});
    end
    function R = getLoadBlock(obj, idx)
      R = obj.loadVec(obj.I{idx}(1):obj.I{idx}(2));
    end
  end
end