classdef Mesh_D
  %% MESH Class
  %
  % General class representing the family of elements and their connectivity
  properties
    element           % Element, reference element
    topology          % MeshTopology
    nodes             % nNxnW, coordinates of vertices
    dimW              % 1x1, dimension of world
    globalSearcher    % GlobalSearcher, searcher of containing element
  end
  methods(Static=true)
    function Mesh(varargin)
      % Constructor
      %
      % Call: m = Mesh(nodes, elem [, dimP])
      %       nodes ... nNxnW, coordinates of vertices
      %       elem  ... nExnV, connectivity of elements
      %       dimP  ... 1x1, space dimension of reference element
      %
      
      %% Tests
      % 1) simple construction 1D
      m = Mesh(linspace(0,1,10)', [(1:9)' (2:10)']);
      if ~isempty(varargin) ,m.show(); pause(); end
      % 2) simple construction 2D
      m = Mesh([0 0; 1 0; 0 1; 1 1], [1 2 3; 4 3 2]);
      if ~isempty(varargin), m.show(); pause(); end
      % 3) simple construction 3D
      m = Mesh([-1 0 0; 1 0 0; 0 1 0; 0 0.5 1], [1 2 3 4]);
      if ~isempty(varargin) ,m.show(); pause(); end
    end
    function evalReferenceMap()
    end
    function evalTrafoInfo()
    end
    function evalInversReferenceMap()
    end
    %
    function evalFunction()
    end
    function integrate()
    end
    %
    function uniformRefine()
    end
    function uniformRefineFast()
    end
    function adaptiveRefine()
    end
    function coarsen()
    end
    %
    function rotate()
    end
    function scale()
    end
    function translate()
    end
    function applyLinearMap()
    end
    %
    function getMeasure()
    end
    function getCenter()
    end
    function getOuterNormal()
    end
    function getOuterNormal2()
    end
    function getDiam()
    end
    function findEntity()
    end
    function findEntityC()
    end
    function isBoundary()
    end
    function getBoundary()
    end
    function isSurface()
    end
    function getSurface()
    end
    function isBoundaryNode()
    end
    function getBoundaryNode()
    end
    %
    function getLocation()
    end
    function showMeshFunction()
    end
    %
    function show()
    end
    %
    function getMesh2D()
    end
    function transformTri2Quad()
    end
    function getTopology()
    end
    function getShapeElement()
    end
  end
end