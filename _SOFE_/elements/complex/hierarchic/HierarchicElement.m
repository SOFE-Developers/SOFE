classdef HierarchicElement < Element
  methods % constructor
    function obj = HierarchicElement(dimension, nV, nB, order)
      obj = obj@Element(dimension, nV, nB, order);
      obj.isLagrange = false;
    end
    function R = getDoFEnum(obj, dim, varargin) % [orient]
      nDoF = obj.doFTuple(dim+1);
      R = (1:nDoF)';
      switch dim
        case 1
          if obj.dimension > 1 && nargin > 2 && varargin{:} < 0
            R(2:2:nDoF) = -R(2:2:nDoF);
          end
        case 2
          if nargin < 3, varargin = {1,1}; end
          if obj.dimension > 2
            if obj.isSimplex()
              R = reshapeTop(floor(sqrt(2*nDoF)):-1:1, R);
              if varargin{2}<0, R = R'; end
              R = R(R>0);
              R = 3*(R-1) + varargin{1};
            else
              n = sqrt(nDoF);
              order1 = repmat((2:n+1).', 1, n); order2 = order1.'; % nxn
              if varargin{3}>0
                sign = varargin{1}.^order1 .* varargin{2}.^order2;
                R = sign(:).*R;
              else
                sign = varargin{1}.^order2 .* varargin{2}.^order1;
                pVec = reshape(1:nDoF,n,n).';
                R = sign(:).*R(pVec(:));
              end
            end
          end
      end
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

