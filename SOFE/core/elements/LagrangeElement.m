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
          if obj.dimension > 1 && ~isempty(varargin) && varargin{1} < 0
            R = (-1)^strcmp(obj.conformity, 'HRot')*R(:,end:-1:1);
          end
        case 2
          if obj.dimension > 2 && nargin > 2
            if obj.isSimplex()
              if nDoF(1)>0
                nCol = floor(sqrt(2*nDoF(1))):-1:1;
                pVec = reshapeTop(nCol, 1:nDoF(1));
                switch varargin{1}
                  case 1
                    % do noting
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
              pVec = 1:nDoF(1);
              if any(cell2mat(varargin(:))<0), pVec = reshape(1:nDoF(1), [], sqrt(nDoF(1))); end
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
  methods % show
    function showLagrangePoints(obj)
      P = obj.getLagrangePoints(obj.dimension, obj.order);
      switch obj.dimension
        case 1
          plot(P,zeros(size(P)),'*');
          hold on
          plot([0 1],[0 0],'r');
          hold off
        case 2
          plot(P(:,1),P(:,2),'*'); axis equal
          hold on
          if obj.isSimplex()
            plot([0 1 0 0]',[0 0 1 0],'r');
          else
            plot([0 1 1 0 0]',[0 0 1 1 0],'r');
          end
          hold off 
        case 3
          fprintf('TODO\n');
      end
    end
  end
end