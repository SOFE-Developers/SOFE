classdef Mesh_D
  %% MESH Class
  %
  % Class representing the mesh, a collections of entities,
  % and their parameterization.
  %
  % SOFE Toolbox.
  % Copyright 2017, Dr. Lars Ludwig
  %
  properties
    element           % <Element, reference element
    topology          % <MeshTopology
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
      % 1) 1D
      m = Mesh(linspace(0,1,10)', [(1:9)' (2:10)']);
      if ~isempty(varargin) ,m.show(); pause(); end
      %
      curve = @(x)[sin(pi*x) cos(pi*x)];
      m2 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      if ~isempty(varargin) ,m2.show(); pause(); end
      %
      curve = @(x)[sin(pi*x) cos(pi*x) x.^2];
      m2 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      if ~isempty(varargin) ,m2.show(); pause(); end
      % 2) 2D
      m = Mesh([0 0; 1 0; 0 1; 1 1], [1 2 3; 4 3 2]);
      if ~isempty(varargin), m.show(); pause(); end
      %
      m = RegularMesh([10 20], [0 1; 0 2]);
      area = @(x)[sin(pi*x(:,1)) cos(pi*x(:,1)) x(:,2)];
      m = Mesh(area(m.nodes), m.topology.getEntity('0'), 2);
      if ~isempty(varargin), m.show(); pause(); end
      % 3) 3D
      m = Mesh([-1 0 0; 1 0 0; 0 1 0; 0 0.5 1], [1 2 3 4]);
      if ~isempty(varargin) ,m.show(); pause(); end
      fprintf('Test ''Mesh.Mesh()'' CHECK\n');
    end
    function evalReferenceMap(varargin)
      % Evaluates function values and first order derivatives of family of 
      % maps from the reference entity in parameter space onto the
      % collection of world entities in physical space.
      %
      % Call: R = m.evalReferenceMap(points, order[, I])
      %       points ... nPxnD, coordinates in parameter space
      %       order  ... 0/1, derivative flag
      %       I      ... nEx1, logical vector selecting a subset of entities
      %       R ... nExnPxnW[xnD], reference map for all points & (sub)entities
      
      %% Tests
      % 1) 1D-->1D/2D/3D
      m = RegularMesh(10, [0 1]);
      R = m.evalReferenceMap(rand(5,1),0);
      assert(norm(size(R)-[10,5])==0);
      R = m.evalReferenceMap(rand(5,1),1);
      assert(norm(size(R)-[10,5])==0);
      %
      curve = @(x)[sin(pi*x) cos(pi*x)];
      m2 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      R = m2.evalReferenceMap(rand(5,1),0);
      assert(norm(size(R)-[10,5,2])==0);
      R = m2.evalReferenceMap(rand(5,1),1);
      assert(norm(size(R)-[10,5,2])==0);
      %
      curve = @(x)[sin(pi*x) cos(pi*x) x.^2];
      m3 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      R = m3.evalReferenceMap(rand(5,1),0);
      assert(norm(size(R)-[10,5,3])==0);
      R = m3.evalReferenceMap(rand(5,1),1);
      assert(norm(size(R)-[10,5,3])==0);
      % 2) 2D-->2D/3D
      m = RegularMesh([10 20], [0 1; 0 2]);
      R = m.evalReferenceMap(rand(5,2),0);
      assert(norm(size(R)-[200,5,2])==0);
      R = m.evalReferenceMap(rand(5,2),1);
      assert(norm(size(R)-[200,5,2,2])==0);
      %
      area = @(x)[sin(pi*x(:,1)) cos(pi*x(:,1)) x(:,2)];
      m3 = Mesh(area(m.nodes), m.topology.getEntity('0'), 2);
      R = m3.evalReferenceMap(rand(5,2),0);
      assert(norm(size(R)-[200,5,3])==0);
      R = m3.evalReferenceMap(rand(5,2),1);
      assert(norm(size(R)-[200,5,3,2])==0);
      % 3) 3D-->3D
      m = RegularMesh([10 5 2], [0 1; 0 2; 0 3]);
      R = m.evalReferenceMap(rand(5,3),0);
      assert(norm(size(R)-[100,5,3])==0);
      R = m.evalReferenceMap(rand(5,3),1);
      assert(norm(size(R)-[100,5,3,3])==0);
      %
      fprintf('Test ''Mesh.evalReferenceMap()'' CHECK\n');
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
  methods(Static=true)
    function testAll()
      fprintf('##########################\n');
      fprintf('Start test suite ''Mesh'':\n');
      fprintf('##########################\n');
      Mesh_D.Mesh();
      Mesh_D.evalReferenceMap();
      fprintf('################################\n');
      fprintf('Test suite ''Mesh'' successful!\n');
      fprintf('################################\n');
    end
  end
end