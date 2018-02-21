classdef EigenSolver < Solver
  properties
    mass
    eigenVal
    eigenVec
  end
  methods % constructor
    function obj = EigenSolver(pde)
      obj = obj@Solver(pde);
      obj.mass = Mass(pde.list{1}.fesTest);
    end
  end
  methods % solve
    function R = compute(obj)
      t = tic; obj.output('Begin assemble ...', 1);
      obj.pde.assemble();
      obj.mass.assemble();
      [freeI, freeJ] = obj.pde.getFreeDoFs();
      fprintf('%d DoFs\n', sum(freeJ));
      %
      obj.output(['... assembled (',num2str(toc(t)),' sec)'], 1);    
      t = tic; obj.output('Begin solve ...', 1);
      A = obj.pde.stiffMat(freeI, freeJ);
      M = obj.mass.stiffMat(freeI, freeJ);
      obj.eigenVec = zeros(size(freeI,1), 20);
      [obj.eigenVec(freeJ,:), obj.eigenVal] = eigs(A,M,20,0,struct('v0',ones(sum(freeI),1)));
      obj.eigenVal = diag(obj.eigenVal);
      obj.output(['... solved (',num2str(toc(t)),' sec)'], 1);
    end
  end
end
