classdef VariationalIntegrator < TimeStep
  % cG(2) method for "M U_t = L - A*U"
  properties
    nK
    freeI, freeJ
    FREEI, FREEJ
    D,G,d,t
  end
  methods % constructor
    function obj = VariationalIntegrator(M0, pde, type, solver)
      obj = obj@TimeStep([], M0, pde, solver);
      obj.nS = 1;
      switch type
        case 'DG1'
          obj.nK = 2;
          obj.D = [3 1;-9 5]/4;
          obj.G = [0 1/2 0; 0 0 1/2];
          obj.d = [1;-1];
          obj.t = [0 1/3 1];
        case 'CG2'
          obj.nK = 2;
          obj.D = [1 1/4;-4 2];
          obj.G = [1/4 1/2 0;-1/2 0 1/2];
          obj.d = [5/4;-2];
          obj.t = [0 1/2 1];
      end
      [obj.freeI, obj.freeJ] = obj.pde.getFreeDoFs();
      obj.FREEI = kron(ones(obj.nK,1), obj.freeI)>0;
      obj.FREEJ = kron(ones(obj.nK,1), obj.freeJ)>0;
    end
  end
  methods % computation
    function R = compute(obj, I, u0)
      A = cell(obj.nK+1,1); b = zeros(numel(u0),obj.nK+1);
      R = cell(obj.nK+1,1);
      GA = cell(1,obj.nK+1);
      obj.pde.assemble();
      A{1} = obj.pde.stiffMat; b(:,1) = obj.pde.loadVec;
      for k = 1:obj.nK
        obj.pde.setState(I(1)+obj.t(k+1)*diff(I), u0);
        obj.pde.assemble();
        A{k+1} = obj.pde.stiffMat; b(:,k+1) = obj.pde.loadVec;
        GA{k+1} = kron(obj.G(:,k+1), A{k+1});
        R{k+1} = obj.pde.getShift();
      end
      lhs = kron(obj.D, obj.M0.stiffMat) + diff(I)*spcell2mat(GA(2:end));
      b = permute(sum(bsxfun(@times,obj.G,permute(b,[3 2 1])),2),[3 1 2]); % nDoF x nB
      rhs = kron(obj.d, obj.M0.stiffMat*u0) + diff(I)*(b(:) - kron(obj.G(:,1), A{1}*u0));
      R = obj.solver.solve(lhs, rhs, obj.FREEI, obj.FREEJ, cell2mat(R));
      R = R(0.5*numel(R)+1:end);
    end
  end
end