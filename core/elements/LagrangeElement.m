classdef LagrangeElement < Element
  properties
    source
    coeffMatrix
  end
  methods % constructor & DoFEnum
    function obj = LagrangeElement(source)
      obj = obj@Element(source.dimension, source.nV, source.nB, source.order);
      obj.isLagrange = true;
      obj.source = source;
      for d = obj.dimension:-1:1
        R = obj.evalFunctionals(d);
        obj.coeffMatrix{d} = eye(size(R,2))/R; % nBxnBs        
        obj.coeffMatrix{d}(abs(obj.coeffMatrix{d})<1e-12) = 0;
      end
    end
    function R = getDoFEnum(obj, dim, varargin) % [orient]
      nDoF = obj.doFTuple(:,dim+1);
      R = reshape((1:prod(nDoF))', [], nDoF(1));
      switch dim
        case 1
          if obj.dimension > 1 && nargin > 2 && varargin{1} < 0
            R = (-1)^strcmp(obj.conformity, 'HRot')*R(:,end:-1:1);
          end
        case 2
          if obj.dimension > 2 && nargin > 2
            if obj.isSimplex()
              if nDoF(1)>0
                nCol = floor(sqrt(2*nDoF(1))):-1:1;
                pVec = reshapeTop(nCol, 1:nDoF(1));
                switch varargin{1}
                  case 2
                    pVec = rot90(pVec);
                    pVec = rot90(reshapeTop(nCol,pVec(pVec>0)));
                  case 3
                    pVec = rot90(pVec);
                end
                if varargin{2} < 0, pVec = reshapeTop(nCol,pVec(pVec>0))'; end                
                R = R(:,pVec(pVec>0));
              end
            else
              pVec = reshape(1:nDoF(1), [], sqrt(nDoF(1)));
              if varargin{1}<0, pVec = pVec(end:-1:1,:); end
              if varargin{2}<0, pVec = pVec(:,end:-1:1); end
              if varargin{3}<0, pVec = pVec.'; end
              R = R(:,pVec);
            end
          end
      end
      R = R(:);
    end
  end
  methods % evaluate
    function B = evalBasis(obj, points, order)
      basis = obj.source.evalBasis(points, order); % nBsxnPxnCx[...]
      sizeVec = size(basis); % nBsxnPxnCx[...]
      sizeVec(1) = size(obj.coeffMatrix{size(points,2)},1); % nBxnPxnCx[...]
      B = reshape(obj.coeffMatrix{size(points,2)}*reshape(basis, size(basis,1), []), sizeVec);
    end
  end
end