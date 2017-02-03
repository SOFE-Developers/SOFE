classdef PpL < LagrangeElement
  methods % constructor
    function obj = PpL(dim, order)
      obj = obj@LagrangeElement(Pp(dim,order));
      if order == 0
        obj.doFTuple = zeros(1,dim+1);
        obj.doFTuple(dim+1) = 1;
      else
        for i = 0:dim
          obj.doFTuple(i+1) = prod(order - (1:i))/factorial(i);
        end
      end
      obj.conformity = 'H1';
    end
    function R = evalFunctionals(obj, dim)
      p = obj.order;
      p1d = linspace(0,1,p+1)';
      p1d = p1d(2:p);
      switch dim
        case 1
          if obj.order == 0
            points = 0.5;
          else
            points = [0; 1];
            points = [points; p1d];
          end
        case 2
          points = [0 0; 1 0; 0 1];
          points = [points; [p1d zeros(p-1,1)]; [1-p1d p1d]; [zeros(p-1,1) p1d]];
          [py,px] = meshgrid(p1d,p1d);
          [iy,ix] = meshgrid(1:p-1, 1:p-1);
          I = ix+iy < p;
          points = [points; [px(I) py(I)]];
        case 3
          zz = zeros(p-1,1);
          points = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
          % edges
          points = [points; ...
                    [p1d zz zz]; [1-p1d p1d zz]; [zz p1d zz]; ...
                    [zz zz p1d]; [1-p1d zz p1d]; [zz 1-p1d p1d]];
          % faces
          [py,px] = meshgrid(p1d,p1d);
          [iy,ix] = meshgrid(1:p-1, 1:p-1);
          I = ix+iy < p;
          zz = zeros(sum(I(:)),1);
          points = [points; ...
                    [px(I) py(I) zz]; [px(I) zz py(I)]; ...
                    [1-px(I)-py(I) px(I) py(I)]; [zz px(I) py(I)]];
          % inner
          [pz,py,px] = meshgrid(p1d,p1d,p1d);
          [iz,iy,ix] = meshgrid(1:p-1, 1:p-1,1:p-1);
          I = ix+iy+iz < p;
          points = [points; [px(I) py(I) pz(I)]];
      end
      R = obj.source.evalBasis(points, 0); % nBxnP
    end
  end
end