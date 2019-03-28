classdef LagrangeElement < Element
  properties
    source
    coeffMatrix
  end
  methods % constructor & DoFEnum
    function obj = LagrangeElement(source, varargin) % [lagPLoc]
      obj = obj@Element(source.dimension, source.nV, source.nB, source.order);
      obj.isLagrange = true;
      obj.source = source;
      for d = obj.dimension:-1:1
        R = obj.evalFunctionals(d, varargin{:});
        obj.coeffMatrix{d} = eye(size(R,2))/R; % nBxnBs        
        obj.coeffMatrix{d}(abs(obj.coeffMatrix{d})<1e-12) = 0;
      end
    end
    function R = getDoFEnum(obj, dim, varargin) % [orient]
      nDoF = obj.doFTuple(:,dim+1);
      R = reshape((1:prod(nDoF))', [], nDoF(1));
      switch dim
        case 1
          if obj.dimension > 1 && ~isempty(varargin) && varargin{1}==2
            R = (-1)^strcmp(obj.conformity, 'HRot')*R(:,end:-1:1);
          end
        case 2
          if obj.dimension > 2 && ~isempty(varargin)
            if obj.isSimplex()
              if nDoF(1)>0
                nCol = floor(sqrt(2*nDoF(1))):-1:1;
                pVec = reshapeTop(nCol, 1:nDoF(1));
                switch ceil(0.5*varargin{1})
                  case 1
                    % do noting
                  case 2
                    pVec = rot90(pVec);
                    pVec = rot90(reshapeTop(nCol,pVec(pVec>0)));
                  case 3
                    pVec = rot90(pVec);
                end
                if mod(varargin{1},2)==0, pVec = reshapeTop(nCol,pVec(pVec>0))'; end
                R = R(:,pVec(pVec>0));
              end
            else
              pVec = reshape(1:nDoF(1), [], sqrt(nDoF(1)));
              flags = [1 1 1;1 1 -1; 1 -1 -1; -1 1 1; -1 1 -1; 1 -1 1; -1 -1 1; -1 -1 -1];
              if flags(varargin{1},3)<0
                if flags(varargin{1},1)<0, pVec = pVec(:,end:-1:1); end
                if flags(varargin{1},2)<0, pVec = pVec(end:-1:1,:); end
                pVec = pVec';
              else
                if flags(varargin{1},1)<0, pVec = pVec(end:-1:1,:); end
                if flags(varargin{1},2)<0, pVec = pVec(:,end:-1:1); end
              end
              R = R(:,pVec);
            end
          end
      end
      R = R(:);
    end
  end
  methods % evaluate
    function B = evalBasis(obj, points, order)
      if isempty(points), B = 1; return; end
      basis = obj.source.evalBasis(points, order); % nBsxnPxnCx[...]
      sizeVec = size(basis); % nBsxnPxnCx[...]
      sizeVec(1) = size(obj.coeffMatrix{size(points,2)},1); % nBxnPxnCx[...]
      B = reshape(obj.coeffMatrix{size(points,2)}*reshape(basis, size(basis,1), []), sizeVec);
    end
  end
  methods
    function R = getProlongationCoeff(obj, varargin)
      [P0,D0] = obj.getLagrangeFunctionals(0);
      F0 = sum(bsxfun(@times,obj.evalBasis(P0,0),permute(D0,[3 1 2])),3)';
      R = cell(1,2^obj.dimension);
      for k = 1:numel(R)
        [P1,D1] = obj.getLagrangeFunctionals(k);
        R{k} = F0 \ sum(bsxfun(@times,obj.evalBasis(P1,0),permute(D1,[3 1 2])),3)';
      end
    end
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
            plot([0 1 0 0],[0 0 1 0],'r');
          else
            plot([0 1 1 0 0],[0 0 1 1 0],'r');
          end
          hold off 
        case 3
          plot3(P(:,1),P(:,2),P(:,3),'*'); axis equal
          hold on
          if obj.isSimplex()
            plot3([0 1 0 0 0 1 0 0],[0 0 1 0 0 0 0 1],[0 0 0 0 1 0 1 0],'r');
          else
            plot3([0 1 1 0 0 0 1 1 1 1 1 1 0 0 0 0], ...
                  [0 0 1 1 0 0 0 0 0 1 1 1 1 1 1 0], ...
                  [0 0 0 0 0 1 1 0 1 1 0 1 1 0 1 1],'r');
          end
          hold off 
      end
    end
  end
end