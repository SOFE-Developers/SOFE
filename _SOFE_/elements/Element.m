classdef Element < SOFEClass
  properties
    dimension
    nV, nB
    order
    doFTuple
    conformity
    isLagrange
  end
  methods % constructor
    function obj = Element(dimension, nV, nB, order)
      obj.dimension = dimension;
      obj.nV = nV;
      obj.nB = nB;
      obj.order = order;
    end
    function R = getDoFEnum(obj, dim, varargin) % [orient]
      if size(obj.doFTuple,1) == 1, obj.doFTuple(2,:) = 1; end
      nDoF = obj.doFTuple(:,dim+1);
      R = (1:prod(nDoF))';
    end
  end
  methods % evaluation local
    function R = evalBasis(obj, points, order) %#ok<STOUT,INUSL>
      eval(['R = obj.evalD' num2str(order) 'Basis(points);']);
    end
  end
  methods % get
    function R = isSimplex(obj)
      nD = obj.dimension;
      R = (obj.nV(nD)==3 && nD==2) || (obj.nV(nD)==4 && nD==3);
    end
  end
  methods % display
    function show(obj, dim, i)
      if i > obj.nB(dim)
        error('XXX: There are only %d basis functions in dim %d!\n', obj.nB(dim), dim);
      end
      gridFine = linspace(0,1,300)';
      gridCoarse = linspace(0,1,20)';
      clf
      switch dim
        case 1
          Z = obj.evalBasis(gridFine, 0);
          plot(gridFine,Z(i,:));
        case 2
          [Xf,Yf] = meshgrid(gridFine, gridFine);
          [Xc,Yc] = meshgrid(gridCoarse, gridCoarse);
          Z = obj.evalBasis([Xf(:) Yf(:)], 0); % nBxnPxnC
          if obj.nV(2) == 3
            Z(:,Xf+Yf>1) = NaN;
          end
          surf(Xf,Yf,reshape(Z(i,:),size(Xf))); 
          shading interp; view(0,90)
        case 3
          N = 30;
          [X,Y,Z] = meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));
          if obj.nV(3) == 4
            I = X+Y+Z<=1;
            X = X(I); Y = Y(I); Z = Z(I);
          end
          P = [X(:) Y(:) Z(:)];
          W = obj.evalBasis(P, 0);
          plot3k(P, 'ColorData', W(i,:)');
      end
    end
    function quiver2D(obj, i, varargin)
      if i > obj.nB(2)
        error('XXX: There are only %d basis functions in dim %d!\n', obj.nB(2), 2);
      end
      LS = linspace(0,1,20)';
      [X,Y] = meshgrid(LS, LS);
      Z = obj.evalBasis([X(:) Y(:)], 0); % nBxnPxnC
      if obj.nV(2) == 3
        Z(:,X+Y>1) = NaN;
      end
      quiver(X,Y,reshape(Z(i,:,1),size(X)),reshape(Z(i,:,2),size(X)), 2, varargin{:});
      axis([0 1 0 1]); axis equal
    end
    function quiver3D(obj, i, varargin)
      if i > obj.nB(3)
        error('XXX: There are only %d basis functions in dim %d!\n', obj.nB(3), 3);
      end
      LS = linspace(0,1,10)';
      [X,Y,Z] = meshgrid(LS, LS, LS);
      W = obj.evalBasis([X(:) Y(:) Z(:)], 0); % nBxnPxnC
      if obj.nV(2) == 3
        W(:,X+Y+Z>1) = NaN;
        plot3([0 1 0 0 0 1 0 0],[0 0 1 0 0 0 0 1],[0 0 0 0 1 0 1 0]);
        hold on
      end
      quiver3(X,Y,Z,reshape(W(i,:,1),size(X)), ...
                   reshape(W(i,:,2),size(X)), ...
                   reshape(W(i,:,3),size(X)), varargin{:});
      hold off
      axis([0 1 0 1 0 1]); axis equal
    end
  end
  methods(Static = true)
    function dkN = getLegendreFunctions(x, N, k)
      dkN = zeros(numel(x),numel(N));
      for i = 1:numel(N)
        n = N(i);
        if(n<k)
          dkN(:,i) = zeros(size(x));
        else
          dkN(:,i) = prod(n+(1:k))*Element.getJacobi(x,n-k,k,k)/2^k;
        end
      end
    end
    function jac = getJacobi(x, N, alpha, beta)
       jac = zeros(numel(x),numel(N));
       for i = 1:numel(N)
         n = N(i);
         if (n < 0)
            jac(:,i) = zeros(size(x));
         elseif (n == 0)
           jac(:,i) = ones(size(x));
         elseif (n == 1)
           jac(:,i) = 0.5*(alpha - beta + (alpha + beta + 2)*x);
         else
           sumAB = alpha + beta;
           jacPrevPrev = ones(size(x));
           jacPrev = 0.5*(alpha - beta + (alpha + beta + 2)*x);
           for k = 2:n
             a1 =  2*k*(k + sumAB)*(2*k + sumAB - 2);
             a2 = (2*k + sumAB - 1)*(alpha*alpha - beta*beta);
             a3 = (2*k + sumAB - 2)*(2*k + sumAB - 1)*(2*k + sumAB);
             a4 =  2*(k + alpha - 1)*(k + beta - 1)*(2*k + sumAB);
             %
             jac(:,i)   = ((a2 + a3*x).*jacPrev - a4*jacPrevPrev)/a1;
             jacPrevPrev = jacPrev;
             jacPrev = jac(:,i);
           end
         end
       end
    end
    function dkN = getShapeFunctions(x, N, k)
      dkN = zeros(numel(x),numel(N));
      for i = 1:numel(N)
        n = N(i);
        if(n==0 && k==0)
          dkN(:,i) = 0.5*(1-x);
        end
        if(n==0 && k==1)
          dkN(:,i) = -0.5.*ones(size(x));
        end
        if(n==1 && k==0)
          dkN(:,i) = 0.5*(1+x);
        end
        if(n==1 && k==1)
          dkN(:,i) = 0.5.*ones(size(x));
        end
        if(n<k && n>=2)
          dkN(:,i) = zeros(size(x));
        end
        if(n>=k && n>=2)
          dkN(:,i) = (prod(n+(1:k))*Element.getJacobi(x,n-k,k,k) - ...
          prod(n+(1:k)-2)*Element.getJacobi(x,n-k-2,k,k))/(2^k*sqrt(4*n-2));
        end
      end
    end
  end
end

