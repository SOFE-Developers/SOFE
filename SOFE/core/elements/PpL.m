classdef PpL < LagrangeElement
  properties
    lagPLoc
  end
  methods % constructor
    function obj = PpL(dim, order, varargin) % [lagPLoc]
      obj = obj@LagrangeElement(Pp(dim,order), varargin{:});
      obj.doFTuple = zeros(2,dim+1);
      if order == 0
        obj.doFTuple(:,dim+1) = 1;
      else
        for i = 0:dim
          obj.doFTuple(1,i+1) = prod(order - (1:i))/factorial(i);
          obj.doFTuple(2,i+1) = 1;
        end
      end
      obj.conformity = 'H1';
    end
    function R = evalFunctionals(obj, dim, varargin) % [lagPLoc]
      points = obj.getLagrangePoints(dim, obj.order, varargin{:});
      R = obj.source.evalBasis(points, 0); % nBxnP
    end
  end
  methods(Static=true)
    function [R,D] = getLagrangePoints(dim, p, varargin) % [lagPLoc]
      if ~isempty(varargin)
        p1d = varargin{1}(:);
      else
        p1d = linspace(0,1,p+1)';
%         p1d = (1+QuadRule.evalWeightedGaussPoints(p+1, @(x)1+0*x(:,1),'Lobatto'))/2;
      end
      switch dim
        case 1
          if p == 0
            R = 0.5;
          else
            R = p1d;
          end
        case 2
          p1d = p1d(2:p);
          % edge
          R = [0 0; 1 0; 0 1];
          R = [R; [p1d zeros(p-1,1)]; [1-p1d p1d]; [zeros(p-1,1) p1d]];
          % inner
          [py,px] = meshgrid(p1d,p1d);
          [iy,ix] = meshgrid(1:p-1, 1:p-1);
          I = ix+iy < p;
          R = [R; [px(I) py(I)]];
        case 3
          p1d = p1d(2:p);
          zz = zeros(p-1,1);
          R = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
          % edge
          R = [R; [p1d zz zz]; [1-p1d p1d zz]; [zz p1d zz]; ...
                  [zz zz p1d]; [1-p1d zz p1d]; [zz 1-p1d p1d]];
          % face
          [py,px] = meshgrid(p1d,p1d);
          [iy,ix] = meshgrid(1:p-1, 1:p-1);
          I = ix+iy < p;
          zz = zeros(sum(I(:)),1);
          R = [R; [px(I) py(I) zz]; [px(I) zz py(I)]; ...
                  [1-px(I)-py(I) px(I) py(I)]; [zz px(I) py(I)]];
          % inner
          [pz,py,px] = meshgrid(p1d,p1d,p1d);
          [iz,iy,ix] = meshgrid(1:p-1, 1:p-1,1:p-1);
          I = ix+iy+iz < p;
          R = [R; [px(I) py(I) pz(I)]];
      end
      D = ones(size(R,1),1);
    end
  end
end