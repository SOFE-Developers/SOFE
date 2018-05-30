classdef RungeKuttaMethod < TimeStep
  properties
    nK
    freeI, freeJ
    FREEI, FREEJ
    solverMass
  end
  methods % constructor
    function obj = RungeKuttaMethod(data, M0, pde, solver)
      obj = obj@TimeStep(data, M0, pde, solver);
      obj.nS = 1; obj.nK = size(obj.data.A,1);
      [obj.freeI, obj.freeJ] = obj.pde.getFreeDoFs();
      obj.FREEI = kron(ones(obj.nK,1), obj.freeI)>0;
      obj.FREEJ = kron(ones(obj.nK,1), obj.freeJ)>0;
      obj.solverMass = DirectSol(1);
      obj.solverMass.pde = obj.pde;
    end
  end
  methods % computation
    function R = compute(obj, I, u0)
      A = cell(1,obj.nK); b = cell(1,obj.nK); shift = cell(obj.nK,1);
      for k = 1:obj.nK
        obj.pde.setState(I(1)+obj.data.c(k)*diff(I), u0);
        obj.pde.assemble();
        b{k} = obj.pde.loadVec;
        A{k} = obj.pde.stiffMat;
        shift{k} = obj.pde.getShift();
      end
      dtA = mat2cell(diff(I)*obj.data.A, obj.nK, ones(obj.nK,1));
      L = cellfun(@(P,Q)kron(P, Q), dtA, b, 'UniformOutput', false);
      isExplicit = norm(reshape(triu(obj.data.A),[],1))==0;
      if isExplicit
        L = mat2cell(sum(cell2mat(L),2), numel(shift{1})*ones(obj.nK,1),1);
        M0u0 = obj.M0.stiffMat*u0;
        R = cell(1,obj.nK);
        for s = 1:obj.nK
          rhs = M0u0 + L{s};
          if s>1
            for k = 1:s-1
              if abs(obj.data.A(s,k))>0
                rhs = rhs - (diff(I)*obj.data.A(s,k))*(A{k}*R{k});
              end
            end
          end
          R{s} = obj.solverMass.solve(obj.M0.stiffMat, rhs, obj.freeI, obj.freeJ, shift{s});
        end
      else
        L = kron(ones(obj.nK,1), obj.M0.stiffMat*u0) + sum(cell2mat(L),2);
        S = cellfun(@(P,Q)kron(P, Q), dtA, A, 'UniformOutput', false);
        S = kron(eye(obj.nK), obj.M0.stiffMat) + cell2mat(S);
        R = obj.solver.solve(S, L, obj.FREEI, obj.FREEJ, cell2mat(shift));
        R = mat2cell(reshape(R,[],obj.nK),size(b{1},1),ones(obj.nK,1));
      end
      R = cellfun(@(a1,a2,a3)a1 - a2*a3, b,A,R, 'UniformOutput', false);
      R = obj.M0.stiffMat*u0 + diff(I)*(cell2mat(R)*obj.data.b');
      obj.pde.setState(I(2), u0);
      R = obj.solverMass.solve(obj.M0.stiffMat, R, obj.freeI, obj.freeJ, obj.pde.getShift());
    end
  end
end