classdef QuadRule < SOFE
  properties
    order % exact for polynomials of degree ORDER
    points
    weights
  end
  methods % constructor
    function obj = QuadRule(n)
      obj.order = n;
      obj.initData();
    end
  end
  methods
    function [P, W] = getGaussPoints(obj, varargin)
      if nargin > 1, type = varargin{1}; else type = 'Gauss'; end
      [P, W] = QuadRule.evalWeightedGaussPoints(ceil((obj.order + 1)/2), @(x)1+0*x(:,1), type);
    end
  end
  methods (Static = true)
    function [P, w] = evalWeightedGaussPoints(n, W, type)
      % computes abcisses P[1..n] and weights w[1..n] of the n-point Gaussian 
      % quadrature formula for the weight function W(x) with
      % type = 'Gauss'   for Gauss-Legendre quadrature
      %      = 'RadauL'  for Gauss-RadauL quadrature
      %      = 'RadauR'  for Gauss-RadauR quadrature
      %      = 'Lobatto' for Gauss-Lobatto quadrature
      %
      a = zeros(n,1);
      b = zeros(n,1);
      oldc = 1;
      for i=1:n
          c = quadgk(@(x)W(x).*QuadRule.evalOrthogonalPolynomial(i-1,x,a,b).^2,-1,1);
          b(i) = c/oldc;
          a(i) = quadgk(@(x)W(x).*x.*QuadRule.evalOrthogonalPolynomial(i-1,x,a,b).^2,-1,1)/c;
          oldc = c;
      end  
      switch type
          case 'Gauss'
          case 'RadauR'
              [p1 p2] = QuadRule.evalOrthogonalPolynomial(n-1,1,a,b);
              a(n) = 1-b(n)*p2/p1;
          case 'RadauL'
              [p1 p2] = QuadRule.evalOrthogonalPolynomial(n-1,-1,a,b);
              a(n) = -1-b(n)*p2/p1;
          case 'Lobatto'
              [p11 p21] = QuadRule.evalOrthogonalPolynomial(n-1,-1,a,b);
              [p12 p22] = QuadRule.evalOrthogonalPolynomial(n-1, 1,a,b);
              A = [p11 p21; p12 p22];
              l = [-p11;p12];
              res  = A\l;
              a(n) = res(1);
              b(n) = res(2);
      end
      z     = diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);
      [V D] = eig(z);
      P     = diag(D);
      w     = b(1)*V(1,:)'.^2;    
    end
    function [returnPj returnPjm1] = evalOrthogonalPolynomial(j,x,a,b)
        if j==0
            returnPj   = ones(size(x));
            returnPjm1 = zeros(size(x));
        else
            returnPjm1 = 0;
            pj         = 1;
            for i=1:j
                returnPj   = (x-a(i)).*pj-b(i)*returnPjm1;
                returnPjm1 = pj;
                pj         = returnPj;
            end
        end
    end
  end
end

