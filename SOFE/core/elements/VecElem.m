classdef VecElem < Element
  properties
    scalarElement
    components
  end
  methods % constructor
    function obj = VecElem(scalE, varargin) % [components]
      if isempty(varargin), components = scalE.dimension; else, components = varargin{1}; end
      obj = obj@Element(scalE.dimension, scalE.nV, scalE.nB*components, scalE.order);
      obj.components = components;
      obj.scalarElement = scalE;
      obj.doFTuple = [scalE.doFTuple(1,:); components*ones(size(scalE.doFTuple(1,:)))];
      obj.isLagrange = scalE.isLagrange;
      obj.conformity = 'H1';
    end
    function R = getDoFEnum(obj, dim, varargin) % [orient]      
      tmp = obj.scalarElement.getDoFEnum(dim, varargin{:});
      R = zeros(obj.components*size(tmp,1),1);
      for d = 1:obj.components
        R(d:obj.components:end) = sign(tmp).*(obj.components*(abs(tmp)-1) + d);
      end
    end
  end
  methods % evaluation
    function B = evalD0Basis(obj, points)
      nC = obj.components;
      Bscalar = obj.scalarElement.evalBasis(points, 0);
      sz = size(Bscalar);
      B = zeros(nC*sz(1), sz(2), nC);
      for d = 1:nC
        B(d:nC:end,:,d) = Bscalar;
      end
    end
    function B = evalD1Basis(obj, points)
      nD = size(points, 2);
      nC = obj.components;
      Bscalar = obj.scalarElement.evalBasis(points, 1);
      B = zeros(obj.nB(nD), size(points,1), nC, nD);
      for d = 1:nC
        B(d:nC:end,:,d,:) = Bscalar;
      end
    end
  end
  methods % Lagrange points
    function [R,D] = getLagrangePoints(obj, dim, r)
      [R,D] = obj.scalarElement.getLagrangePoints(dim, r); % nPxnD
      R = kron(R,ones(obj.components,1));
      D = kron(D,eye(obj.dimension));
    end
  end
end