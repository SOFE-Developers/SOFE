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
    function Mesh()
      % Mesh class constructor
      %
      % Call:
      % m = Mesh(nodes, elem [, dimP])
      %
      % nodes ... nNxnW, coordinates of vertices
      % elem  ... nExnV, connectivity of elements
      % dimP  ... 1x1, space dimension of reference element
      %
    end
    function evalReferenceMap()
      % Evaluates function values and first order derivatives of family of 
      % maps from the reference entity in parameter space onto the
      % collection of world entities in physical space.
      %
      % Call:
      % R = m.evalReferenceMap(x, order[, I])
      %
      % x ... nPxnD, coordinates in parameter space OR
      % x ... {nPxnD, nP}, pairs of coordinates in parameter space and hosting element
      % order  ... 0/1, derivative flag
      % I      ... nEx1, logical vector selecting a subset of entities
      % R      ... nExnPxnW[xnD], reference map for all points & (sub)entities
      %
    end
    function evalTrafoInfo()
      % Evaluates reference maps and derives necessary quanitities for
      % transformation of integrals.
      %
      % Call: 
      % [R, invR, jacR] = m.evalTrafoInfo(x[, I])
      %
      % x ... nPxnD, coordinates in parameter space OR
      % x ... {nPxnD, nP}, pairs of coordinates in parameter space and hosting element
      % I      ... nEx1, logical vector selecting a subset of entities
      % R      ... nExnPxnWxnD, gradient of reference map for all points & (sub)entities
      % invR   ... nExnPxnW[xnD], inverse of gradient
      % jacR   ... nExnPxnW[xnD], determinant of gradient
      %
    end
    function evalInversReferenceMap()
      % Finds hosting cell K and inverse reference map related to K for list of
      % global points in physical space.
      %
      % Call:
      % R = m.evalInversReferenceMap(x)
      %
      % x ... nPxnW, global coordinates of query points
      % R ... {L,H} WITH
      % L ... nPxnD, local coordinates of query points AND
      % H ... nPx1, index of hosting cell
      %
    end
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
  methods(Static=true) % tests
    function testAll()
      fprintf('##########################\n');
      fprintf('Start test suite ''Mesh'':\n');
      fprintf('##########################\n');
      Mesh_D.Mesh_T();
      Mesh_D.evalReferenceMap_T();
      Mesh_D.evalTrafoInfo_T();
      fprintf('################################\n');
      fprintf('Test suite ''Mesh'' successful!\n');
      fprintf('################################\n');
    end
    function Mesh_T(varargin)
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
    function evalReferenceMap_T(varargin)
      dx = 1e-8;
      % 1) 1D-->1D/2D/3D
      x = 0.5*rand(5,1);
      m = RegularMesh(10, [0 1]);
      R = m.evalReferenceMap(x,0);
      assert(norm(size(R)-[10,5])==0);
      dR = m.evalReferenceMap(x,1);
      assert(norm(size(dR)-[10,5])==0);
      small = (m.evalReferenceMap(x+dx,0) - R)/dx - dR;
      assert(norm(small(:))<1e-6);
      %
      curve = @(x)[sin(pi*x) cos(pi*x)];
      m2 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      R = m2.evalReferenceMap(x,0);
      assert(norm(size(R)-[10,5,2])==0);
      dR = m2.evalReferenceMap(x,1);
      assert(norm(size(dR)-[10,5,2])==0);
      small = (m2.evalReferenceMap(x+dx,0) - R)/dx - dR;
      assert(norm(small(:))<1e-6);
      %
      curve = @(x)[sin(pi*x) cos(pi*x) x.^2];
      m3 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      R = m3.evalReferenceMap(x,0);
      assert(norm(size(R)-[10,5,3])==0);
      dR = m3.evalReferenceMap(x,1);
      assert(norm(size(dR)-[10,5,3])==0);
      small = (m3.evalReferenceMap(x+dx,0) - R)/dx - dR;
      assert(norm(small(:))<1e-6);
      % 2) 2D-->2D/3D
      x = 0.5*rand(5,2);
      m = RegularMesh([10 20], [0 1; 0 2]);
      R = m.evalReferenceMap(x,0);
      assert(norm(size(R)-[200,5,2])==0);
      dR = m.evalReferenceMap(x,1);
      assert(norm(size(dR)-[200,5,2,2])==0);
      small = (m.evalReferenceMap(x+dx*[1 0],0) - R)/dx - dR(:,:,:,1);
      assert(norm(small(:))<1e-6);
      small = (m.evalReferenceMap(x+dx*[0 1],0) - R)/dx - dR(:,:,:,2);
      assert(norm(small(:))<1e-6);
      %
      area = @(x)[sin(pi*x(:,1)) cos(pi*x(:,1)) x(:,2)];
      m3 = Mesh(area(m.nodes), m.topology.getEntity('0'), 2);
      R = m3.evalReferenceMap(x,0);
      assert(norm(size(R)-[200,5,3])==0);
      dR = m3.evalReferenceMap(x,1);
      assert(norm(size(dR)-[200,5,3,2])==0);
      small = (m3.evalReferenceMap(x+dx*[1 0],0) - R)/dx - dR(:,:,:,1);
      assert(norm(small(:))<1e-6);
      small = (m3.evalReferenceMap(x+dx*[0 1],0) - R)/dx - dR(:,:,:,2);
      assert(norm(small(:))<1e-6);
      % 3) 3D-->3D
      x = 0.5*rand(5,3);
      m = RegularMesh([10 5 2], [0 1; 0 2; 0 3]);
      R = m.evalReferenceMap(x,0);
      assert(norm(size(R)-[100,5,3])==0);
      dR = m.evalReferenceMap(x,1);
      assert(norm(size(dR)-[100,5,3,3])==0);
      small = (m.evalReferenceMap(x+dx*[1 0 0],0) - R)/dx - dR(:,:,:,1);
      assert(norm(small(:))<1e-6);
      small = (m.evalReferenceMap(x+dx*[0 1 0],0) - R)/dx - dR(:,:,:,2);
      assert(norm(small(:))<1e-6);
      small = (m.evalReferenceMap(x+dx*[0 0 1],0) - R)/dx - dR(:,:,:,3);
      assert(norm(small(:))<1e-6);
      %
      fprintf('Test ''Mesh.evalReferenceMap()'' CHECK\n');
    end
    function evalTrafoInfo_T(varargin)
      % 1) 1D-->1D/2D/3D
      x = 0.5*rand(5,1);
      m = RegularMesh(10, [0 1]);
      [R, invR, jacR] = m.evalTrafoInfo(x,1:2:10);
      assert(norm(size(R)-[5,5])==0);
      small = R-1./invR;
      assert(norm(small(:))<1e-14);
      small = R-jacR;
      assert(norm(small(:))<1e-14);
      %
      curve = @(x)[sin(pi*x) cos(pi*x)];
      m2 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      [R, invR, jacR] = m2.evalTrafoInfo(x,1:2:10);
      assert(norm(size(R)-[5,5,2])==0);
      small = permute(R,[1 2 4 3])./sum(R.*R,3)-invR;
      assert(norm(small(:))<1e-14);
      sum(R.*R,3).^0.5 - jacR;
      assert(norm(small(:))<1e-14);
      %
      curve = @(x)[sin(pi*x) cos(pi*x) x.^2];
      m3 = Mesh(curve(m.nodes), m.topology.getEntity('0'), 1);
      [R, invR, jacR] = m3.evalTrafoInfo(x,1:2:10);
      assert(norm(size(R)-[5,5,3])==0);
      small = permute(R,[1 2 4 3])./sum(R.*R,3)-invR;
      assert(norm(small(:))<1e-14);
      sum(R.*R,3).^0.5 - jacR;
      assert(norm(small(:))<1e-14);
      % 2) 2D-->2D/3D
      x = 0.5*rand(5,2);
      m = RegularMesh([10 20], [0 1; 0 2]);
      [R, invR, jacR] = m.evalTrafoInfo(x);
      assert(norm(size(R)-[200,5,2,2])==0);
      small = inv(squeeze(R(1,1,:,:))) - squeeze(invR(1,1,:,:));
      assert(norm(small(:))<1e-14);
      small = det(squeeze(R(1,1,:,:))) - squeeze(jacR(1,1,:,:));
      assert(norm(small(:))<1e-14);
      %
      area = @(x)[sin(pi*x(:,1)) cos(pi*x(:,1)) x(:,2)];
      m3 = Mesh(area(m.nodes), m.topology.getEntity('0'), 2);
      [R, invR, jacR] = m3.evalTrafoInfo(x);
      assert(norm(size(R)-[200,5,3,2])==0);
      r = squeeze(R(1,1,:,:));
      small = (r*inv(r'*r))' - squeeze(invR(1,1,:,:));
      assert(norm(small(:))<1e-14);
      small = det(r'*r).^0.5 - squeeze(jacR(1,1,:,:));
      assert(norm(small(:))<1e-14);
      % 3) 3D-->3D
      x = 0.5*rand(5,3);
      m = RegularMesh([10 5 2], [0 1; 0 2; 0 3]);
      [R, invR, jacR] = m.evalTrafoInfo(x);
      assert(norm(size(R)-[100,5,3,3])==0);
      small = inv(squeeze(R(1,1,:,:))) - squeeze(invR(1,1,:,:));
      assert(norm(small(:))<1e-14);
      small = det(squeeze(R(1,1,:,:))) - squeeze(jacR(1,1,:,:));
      assert(norm(small(:))<1e-14);
      %
      fprintf('Test ''Mesh.evalTrafoInfo()'' CHECK\n');
    end
  end
end