classdef DirectSolver < Solver
  properties
    isCaching = false;
    cache
  end
  methods % constructor
    function obj = DirectSolver(pde, varargin) % [isLU]
      obj = obj@Solver(pde);
      if nargin > 1, obj.isCaching = varargin{1}; end
    end
    function notify(obj)
      obj.cache = [];
    end
  end
  methods % solve
    function R = compute(obj)
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      [freeI, freeJ] = obj.pde.getFreeDoFs();
      fprintf('%d DoFs\n', size(obj.pde.loadVec,1));
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);
      %
      t = tic; obj.output('Begin solve ...', 1);
      if ~obj.pde.createSys 
        error('System must be created');
      end
      R = obj.pde.getShift();
      b = obj.pde.loadVec - obj.pde.stiffMat*R;
      b = b(freeI);
      A = obj.pde.stiffMat(freeI, freeJ);
      %
      if obj.isCaching
        if isempty(obj.cache)
          [obj.cache.L,obj.cache.U,obj.cache.P,obj.cache.Q,obj.cache.R] = lu(A);
        end
        R(freeJ) = obj.cache.Q*(obj.cache.U\(obj.cache.L\(obj.cache.P*(obj.cache.R\b))));
      else
        R(freeJ) = A\b;
      end
      obj.solution = R;
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
    function R = compute2(obj) % extension to hanging nodes
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      fprintf('%d DoFs\n', size(obj.pde.loadVec,1));
      %
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);    
      t = tic; obj.output('Begin solve ...', 1);
      if ~obj.pde.createSys 
        error('System must be created');
      end
      R = obj.pde.getShift();
      b = obj.pde.loadVec - obj.pde.stiffMat*R;
      PP = obj.pde.fesTrial{1}.getProjector();
      [freeI, freeJ] = obj.pde.getFreeDoFs();
      if isempty(PP)
        PP = speye(numel(freeI));
        P = PP; Q = P;
        P(~freeI,:) = []; Q(:,~freeJ) = [];
      else
        assert(norm(size(PP)-size(obj.pde.stiffMat))==0,'size error');
        P = PP; Q = PP';
        I = (sum(P,2)==0);
        P(~freeI | I,:) = [];
        Q(:, ~freeJ | I) = [];
      end
      if obj.isCaching
        if isempty(obj.cache)
          [obj.cache.L,obj.cache.U,obj.cache.P,obj.cache.Q,obj.cache.R] = lu(P*obj.pde.stiffMat*Q);
        end
        Rtmp = obj.cache.Q*(obj.cache.U\(obj.cache.L\(obj.cache.P*(obj.cache.R\(P*b)))));
      else
        Rtmp = (P*obj.pde.stiffMat*Q)\(P*b);
      end
      obj.solution = PP'*R + Q*Rtmp;
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
