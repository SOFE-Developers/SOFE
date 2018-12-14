classdef VariationalIntegrator < TimeStep
  properties
    type
    D,G,L,d,t,pad
    basis0I, basis0J
    basisQI, dBasisQI
    basisQJ, dBasisQJ
    fesTrial, fesTest
  end
  methods % constructor and init
    function obj = VariationalIntegrator(M0, pde, type, solver, varargin) % [weightfunc]
      obj = obj@TimeStep([], M0, pde, solver);
      obj.nS = 1;
      order = str2double(type(3:end));
      obj.type = type(1:2);
      obj.init(order, varargin{:});
    end
    function init(obj, order, varargin) % [weightfunc]
      obj.nK = order+1;
      switch obj.type
        case {'dG','DG'}
          quadType = 'RadauR';
        case {'cG','CG'}
          quadType = 'Lobatto';
        otherwise
          error('unknown type: %s', obj.type);
      end
      if ~isempty(varargin), wf = varargin{1}; else, wf = @(t)1+0*t; end
      [qp,w]=QuadRule.evalWeightedGaussPoints(obj.nK,wf,quadType); qp=(qp+1)/2;w=w/2;
      e = PpL(1, order, qp);
      switch obj.type
        case {'dG','DG'}
          eTest = PpL(1, order, qp);
        case {'cG','CG'}
          eTest=PpL(1, order-1);
      end
      e.quadRule{1}.points=qp; e.quadRule{1}.weights = w;
      eTest.quadRule{1}.points = qp; eTest.quadRule{1}.weights = w;
      obj.t   = e.quadRule{1}.points';
      obj.fesTrial = FESpace(RegularMesh(1,[0,1]), e);
      obj.fesTest  = FESpace(RegularMesh(1,[0,1]), eTest);
      mass  = OpIdId(1, 0, obj.fesTrial,obj.fesTest);
      stiff = OpGradId(1, obj.fesTrial, obj.fesTest);
      mass.assemble();
      stiff.assemble();
      load = FcId(1, obj.fesTest, 0); load.assemble();
      obj.basis0I = obj.fesTest.element.evalBasis(0, 0); % nB
      obj.basis0J = obj.fesTrial.element.evalBasis(0, 0); % nB
      obj.basisQI = obj.fesTest.element.evalBasis(qp,0); % nBxnP
      obj.basisQJ = obj.fesTrial.element.evalBasis(qp,0); % nBxnP
      obj.dBasisQI = obj.fesTest.element.evalBasis(qp,1); % nBxnP
      obj.dBasisQJ = obj.fesTrial.element.evalBasis(qp,1); % nBxnP
      obj.basisQI(abs(obj.basisQI)<1e-14) = 0;
      obj.basisQJ(abs(obj.basisQJ)<1e-14) = 0;
      obj.dBasisQI(abs(obj.dBasisQI)<1e-14) = 0;
      obj.dBasisQJ(abs(obj.dBasisQJ)<1e-14) = 0;
      obj.G = full(mass.matrix); obj.G(abs(obj.G)<1e-14) = 0;
      obj.D = full(stiff.matrix); obj.D(abs(obj.D)<1e-14) = 0;
      obj.L = load.matrix;
      switch obj.type
        case 'dG'
          obj.D   = obj.D + obj.basis0I*obj.basis0J';
          obj.d   = obj.basis0I;
          obj.pad = [];
        case 'cG'
          obj.d   = [];
          obj.pad = obj.basis0J';
      end
    end
  end
  methods % computation
    function R = compute(obj, I, u0) % variable data, constant M0
      A = cell(1,obj.nK); 
      b = zeros(numel(u0),obj.nK);
      shift = cell(obj.nK,1);
      for k = 1:obj.nK
        obj.pde.setState(I(1)+obj.t(k)*diff(I), u0);
        obj.pde.assemble();
        A{k}   = kron(obj.G(:,k), obj.pde.stiffMat); 
        b(:,k) = obj.pde.loadVec;
        shift{k}   = obj.pde.getShift();
      end
      lhs = kron(obj.D, obj.M0.stiffMat) + diff(I)*spcell2mat(A);
      %
      w = diff(I)*obj.fesTest.element.quadRule{1}.weights;
      rhs = bsxfun(@times, b, permute(obj.basisQI,[3 2 1])); % nDoFxnPxnB
      rhs = sum(bsxfun(@times, rhs, w'), 2); % nDoFx1xnB
      rhs = reshape(permute(rhs, [1 3 2]), [], 1); % nDoF*nB
      %
      if numel(obj.d>0) % add jump contribution, dG
        rhs = rhs + kron(obj.d, obj.M0.stiffMat*u0);
      end
      if numel(obj.pad>0) % add continuity condition, cG
        lhs = [lhs; kron(obj.pad, obj.M0.stiffMat)];
        rhs = [rhs; obj.M0.stiffMat*u0];
      end
      R = obj.solver.solve(lhs, rhs, obj.FREEI, obj.FREEJ, cell2mat(shift));
      R = reshape(R,[],obj.nK);
    end
    function R = compute_(obj, I, u0) % constant data
      lhs = kron(obj.D, obj.M0.stiffMat) + diff(I)*kron(obj.G, obj.pde.stiffMat);
      rhs = diff(I)*kron(obj.L, obj.pde.loadVec);
      shift = kron(ones(obj.nK,1),obj.pde.getShift());
      %
      if numel(obj.d>0) % add jump contribution, dG
        rhs = rhs + kron(obj.d, obj.M0.stiffMat*u0);
      end
      if numel(obj.pad>0) % add continuity condition, cG
        lhs = [lhs; kron(obj.pad, obj.M0.stiffMat)];
        rhs = [rhs; obj.M0.stiffMat*u0];
      end
      R = obj.solver.solve(lhs, rhs, obj.FREEI, obj.FREEJ, shift);
      R = reshape(R,[],obj.nK);
    end
    function R = compute__(obj, I, u0) % variable M0
      A = cell(1,obj.nK);
      B = cell(1,obj.nK);
      DD = cell(1,obj.nK);
      M = cell(1,obj.nK+1);
      b = zeros(numel(u0),obj.nK);
      R = zeros(numel(u0),obj.nK);
      w = obj.fesTest.element.quadRule{1}.weights;
      %
      obj.pde.setState(I(1), u0);
      M{1}   = obj.M0.stiffMat;
      for k = 1:obj.nK
        obj.pde.setState(I(1)+obj.t(k)*diff(I), u0);
        obj.M0.setState(I(1)+obj.t(k)*diff(I), u0);
        obj.pde.assemble(); obj.M0.assemble();
        A{k} = obj.pde.stiffMat;
        B{k} = ((w(k)*diff(I))*obj.basisQI(:,k)).*obj.basisQJ(:,k)';
        DD{k} = (w(k)*obj.basisQI(:,k)).*obj.dBasisQJ(:,k)';
        M{k+1} = obj.M0.stiffMat;
        b(:,k) = obj.pde.loadVec;
        R(:,k) = obj.pde.getShift();
      end
      % get shift
      lhs = sparse(size(B{1},1)*size(A{1},1), size(B{1},2)*size(A{1},2));
      for k = 1:obj.nK
        lhs = lhs + kron(B{k}, speye(size(M{k+1},1),size(M{k+1},2)));  
      end
      rhs = bsxfun(@times, R, permute(obj.basisQI,[3 2 1])); % nDoFxnPxnB
      rhs = sum(bsxfun(@times, rhs, w'), 2); % nDoFx1xnB
      rhs = diff(I)*reshape(permute(rhs, [1 3 2]), [], 1); % nDoF*nB
      R = lhs\rhs;
      %
      switch obj.type
        case 'dG'
          lhs = kron(obj.basis0I.*obj.basis0J', M{1});
        case 'cG'
          lhs = sparse(size(B{1},1)*size(A{1},1), size(B{1},2)*size(A{1},2));
      end
      for k = 1:obj.nK
        lhs = lhs + kron(B{k}, A{k}) + kron(DD{k}, M{k+1});  
      end
      rhs = bsxfun(@times, b, permute(obj.basisQI,[3 2 1])); % nDoFxnPxnB
      rhs = bsxfun(@times, rhs, diff(I)*w'); % nDoFx1xnB
      rhs = reshape(permute(sum(rhs,2), [1 3 2]), [], 1); % nDoF*nB
      switch obj.type
        case 'dG'
          rhs = rhs + kron(obj.basis0I, M{1}*u0);
        case 'cG'
          lhs = [lhs; kron(obj.basis0J', M{1})];
          rhs = [rhs; M{1}*u0];
      end
      R = obj.solver.solve(lhs, rhs, obj.FREEI, obj.FREEJ, R(:));
      R = reshape(R,[],obj.nK);
    end
    function R = eval(obj, H, tLoc)
      R = permute(obj.fesTrial.evalDoFVectorLocal(H', tLoc(:), [], 0), [3 2 1]); % NxnP
    end
  end
end