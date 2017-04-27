classdef PDE < SOFE
  properties
    nEq
    lhs, rhs
    %
    stiffMat, loadVec, solution
    shift, fDoFsTrial, fDoFsTest
    %
    mesh
    fesTest,fesTrial
    I,J, nDoF
    %
    time, state
    nonLin = false;
    stateChanged
    %
    solver = DirectSolver([]);
  end
  methods % constructor & more
    function obj = PDE(lhs, rhs, varargin)
      obj.lhs = lhs;
      obj.rhs = rhs;
      obj.nEq = numel(rhs);
      if nargin > 2
        obj.nonLin = varargin{1}.nonLin;
      end
      obj.fesTrial = cell(obj.nEq, 1);
      obj.fesTest = cell(obj.nEq, 1);
      for i = 1:obj.nEq
        for j = 1:obj.nEq
          if ~isempty(obj.lhs{i,j})
            if isempty(obj.fesTrial{j})
              obj.fesTrial{j} = obj.lhs{i,j}{1}.fesTrial;
              obj.fesTrial{j}.register(obj);
            end
            if isempty(obj.fesTest{i})
              obj.fesTest{i} = obj.lhs{i,j}{1}.fesTest;
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
      for i = 1:obj.nEq
        dimTest(i) = obj.fesTest{i}.getNDoF();
        dimTrial(i) = obj.fesTrial{i}.getNDoF();
      end
      offsetTrial = [0 cumsum(dimTrial)]; offsetTest = [0 cumsum(dimTest)];
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
      obj.time = newTime;
      for k = 1:obj.nEq
        if isempty(obj.fesTrial{k}.shift) || nargin(obj.fesTrial{k}.shift) > 1
          obj.shift = [];
        end
        if isempty(obj.fesTrial{k}.fixB) || nargin(obj.fesTrial{k}.fixB) > 1
          obj.fDoFsTrial = [];
        end
        if isempty(obj.fesTest{k}.fixB) || nargin(obj.fesTest{k}.fixB) > 1
          obj.fDoFsTest = [];
        end
      end
      obj.stateChanged = true;
    end
    function setState(obj, varargin)
      obj.stateChanged = true;
      if ~obj.nonLin, return; end
      obj.state = cell(obj.nEq, 1); % {nEq}xnExnP
      for j = 1:obj.nEq
        nBlock = obj.fesTrial{j}.nBlock(1);
        [~, w] = obj.fesTrial{j}.getQuadData(0);
        nD = obj.fesTrial{j}.element.dimension;
        nC = size(obj.fesTrial{j}.element.evalBasis(zeros(1,nD),0),3);
        obj.state{j} = zeros(obj.mesh.topology.getNumber(obj.mesh.topology.dimP), numel(w), nC);
        if nargin < 2, continue; end
        U = varargin{1}(obj.J{j}(1):obj.J{j}(2));
        obj.state{j} = obj.fesTrial{j}.evalDoFVector(U, [], 0, 0);
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
      obj.setState(obj.solution);
    end
    function assemble(obj)
      if obj.stateChanged, 
        obj.stateChanged = false;
      else
        return; 
      end
      % lhs
      obj.stiffMat = sparse(obj.I{obj.nEq}(2), obj.J{obj.nEq}(2));
      for i = 1:obj.nEq
        for j = 1:obj.nEq
          if ~isempty(obj.lhs{i,j})
            for k = 1:numel(obj.lhs{i,j})
              obj.lhs{i,j}{k}.notify(obj.time, obj.state);
              obj.lhs{i,j}{k}.assemble();
              %
              blk = obj.lhs{i,j}{k}.matrix;
              blk = [sparse(obj.I{i}(1)-1, size(blk,2)); blk];
              blk(obj.I{i}(2)+1:obj.I{obj.nEq}(2),:) = 0;
              blk = [sparse(size(blk,1), obj.J{j}(1)-1), blk];
              blk(:,obj.J{j}(2)+1:obj.J{obj.nEq}(2)) = 0;
              %
              obj.stiffMat = obj.stiffMat + blk;
            end
          end
        end
      end
      % rhs
      obj.loadVec = zeros(obj.I{obj.nEq}(2), 1);
      for i = 1:obj.nEq
        if ~isempty(obj.rhs{i})
          for k = 1:numel(obj.rhs{i})
            obj.rhs{i}{k}.notify(obj.time, obj.state);
            obj.rhs{i}{k}.assemble();
            idx = obj.I{i}(1):obj.I{i}(2);
            obj.loadVec(idx,:) = obj.loadVec(idx,:) + obj.rhs{i}{k}.vector;
          end
        end
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
    end
    function solve(obj)
      obj.solution = zeros(size(obj.loadVec));
      for k = 1:size(obj.loadVec,2)
        b = obj.loadVec(:,k);
        b = b - obj.stiffMat*obj.shift; 
        obj.solution(~obj.fDoFsTrial, k) = obj.shift(~obj.fDoFsTrial);
        A = obj.stiffMat(obj.fDoFsTest,obj.fDoFsTrial);
        b = b(obj.fDoFsTest);
        x = obj.solver.solve(A, b);
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
    function R = getStiffnessBlock(obj, idxI, idxJ)
      R = obj.stiffMat(obj.I{idxI},obj.J{idxJ});
    end
    function R = getLoadBlock(obj, idx)
      R = obj.loadVec(obj.I{idx}(1):obj.I{idx}(2));
    end
  end
end