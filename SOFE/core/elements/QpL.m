classdef QpL < LagrangeElement
  methods % constructor
    function obj = QpL(dim, p)
      obj = obj@LagrangeElement(Qp(dim,p));
      obj.order = p;
      if p == 0
        obj.doFTuple = zeros(1,dim+1);
        obj.doFTuple(dim+1) = 1;
      else
        for i = 0:dim
          obj.doFTuple(i+1) = (obj.order-1)^i;
        end
      end
      obj.conformity = 'H1';
    end
    function R = evalFunctionals(obj, dim)
      points = obj.getLagrangePoints(dim, obj.order);
      R = obj.source.evalBasis(points, 0); % nBxnP
    end
  end
  methods(Static=true)
    function R = getLagrangePoints(dim, p)
      p = p(1);
      if p == 0
        p1d = 0.5;
      else
%         p1d = linspace(0,1,p+1)';
        p1d = (1+QuadRule.evalWeightedGaussPoints(p+1, @(x)1+0*x(:,1),'Lobatto'))/2;
      end
      p1d = p1d(2:p);
      switch dim
        case 1
          R = [0; 1; p1d];
        case 2
          zz = zeros(p-1,1); oo = ones(p-1,1);
          R = [0 0; 1 0; 0 1; 1 1];
          % faces
          R = [R; [p1d zz]; [p1d oo]; [zz p1d]; [oo p1d]];
          % inner
          [py,px] = meshgrid(p1d,p1d);
          R = [R; [px(:) py(:)]];
        case 3
          zz = zeros(p-1,1); oo = ones(p-1,1);
          R = [0 0 0; 1 0 0; 0 1 0; 1 1 0; ...
               0 0 1; 1 0 1; 0 1 1; 1 1 1];
          % edges
          R = [R; [p1d zz zz]; [p1d oo zz]; [p1d zz oo]; [p1d oo oo]; ...
                  [zz p1d zz]; [zz p1d oo]; [oo p1d zz]; [oo p1d oo]; ...
                  [zz zz p1d]; [oo zz p1d]; [zz oo p1d]; [oo oo p1d]];
          % faces
          [py,px] = meshgrid(p1d,p1d);
          zz = zeros((p-1)^2,1); oo = ones((p-1)^2,1);
          R = [R; [px(:) py(:) zz]; [px(:) py(:) oo]; ...
                  [zz px(:) py(:)]; [oo px(:) py(:)]; ...
                  [px(:) zz py(:)]; [px(:) oo py(:)]];
          % inner
          [pz,py,px] = meshgrid(p1d,p1d,p1d);
          R = [R; [px(:) py(:) pz(:)]];
      end
    end
  end
end