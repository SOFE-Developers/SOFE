classdef TPElem < Element
  properties
    scalarElement
  end
  methods % constructor
    function obj = TPElem(scalE)
      obj = obj@Element(scalE.dimension, scalE.nV, scalE.nB*scalE.dimension, scalE.order);
      obj.scalarElement = scalE;
      obj.doFTuple = [scalE.doFTuple; scalE.dimension*ones(size(scalE.doFTuple))];
      obj.isLagrange = scalE.isLagrange;
      obj.conformity = 'H1';
    end
    function R = getDoFEnum(obj, dim, varargin) % [orient]      
      tmp = obj.scalarElement.getDoFEnum(dim, varargin{:});
      R = zeros(obj.dimension*size(tmp,1),1);
      for d = 1:obj.dimension
        R(d:obj.dimension:end) = sign(tmp).*(obj.dimension*(abs(tmp)-1) + d);
      end
    end
  end
  methods % evaluation
    function B = evalD0Basis(obj, points)
      nD = size(points, 2);
      nW = obj.dimension;
      Bscalar = obj.scalarElement.evalBasis(points, 0);
      B = zeros(obj.nB(nD), size(points,1), nD);
      for d = 1:nW
        B(d:nW:end,:,d) = Bscalar;
      end
    end
    function B = evalD1Basis(obj, points)
      nD = size(points, 2);
      nW = obj.dimension;
      Bscalar = obj.scalarElement.evalBasis(points, 1);
      B = zeros(obj.nB(nD), size(points,1), nD, nD);
      for d = 1:nW
        B(d:nW:end,:,d,:) = Bscalar;
      end
    end
  end
end