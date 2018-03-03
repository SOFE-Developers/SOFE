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
    % evaluation reference map
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
      % I      ... nEx1, logical vector selecting subset of entities
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
      % I      ... nEx1, logical vector selecting subset of entities
      % R      ... nExnPxnWxnD, gradient of reference map for all points & (sub)entities
      % invR   ... nExnPxnW[xnD], inverse of gradient
      % jacR   ... nExnPxnW[xnD], determinant of gradient
      %
    end
    function evalInversReferenceMap()
      % Finds hosting cell K and evaluates inverse reference map related to K 
      % for all global points in physical space contained in list.
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
      % Evaluates (nonlinear) function in local points.
      %
      % Call:
      % R = m.evalFunction(F, x, S[, I])
      %
      % F = F(x[,u,du]) ... function handle
      % x    ... nPxnD, local coordinates
      % S    ... [] OR struct('U','dU') WITH
      % S.U  ... cell(mxn){nExnPxnC}, evaluated state
      % S.dU ... cell(mxn){nExnPxnCxnD}, gradient of evaluated state
      % I    ... nEx1, logical vector selecting subset of entities
      % R    ... nExnPxnCxnD, array function values
      %
    end
    function integrate()
      % Integrates function on (sub)mesh.
      %
      % Call:
      % R = m.integrate(F, quadRule[, I])
      %
      % F = F(x) ... function handle
      % quadRule ... QuadRule, quadrature rule class
      % I    ... nEx1, logical vector selecting subset of entities
      % R    ... 1x1, value of integral
      %
    end
    % refinement
    function uniformRefine()
      % Refines all cells uniformly in 2^nD subcells (red refinement)
      %
      % Call:
      % m.uniformRefine([N])
      %
      % N ... 1x1, number of refinements, default N = 1
      %
    end
    function uniformRefineFast()
      % Refines all cells uniformly in 2^nD subcells (red refinement)
      % using fast update of connectivity structures
      %
      % Call:
      % m.uniformRefineFast([N])
      %
      % N ... 1x1, number of refinements, default N = 1
      %
    end
    function adaptiveRefine()
      % Refines all cells in given subdomain adaptively by bisection
      %
      % Call:
      % m.adaptiveRefine(loc[,N])
      %
      % loc ... function handle, location of cells to be refined
      % N   ... 1x1, number of refinements, default N = 1
      %
    end
    function coarsen()
      % Coarsens all cells in given subdomain that previously have been refined
      % adaptively.
      %
      % Call:
      % m.coarsen(loc[,N])
      %
      % loc ... function handle, location of cells to be refined
      % N   ... 1x1, number of refinements, default N = 1
      %
    end
    % mesh shape manipulation
    function rotate()
    end
    function scale()
    end
    function translate()
    end
    function applyLinearMap()
    end
    % mesh information
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
    % mesh visualization
    function show()
    end
    % mesh generation
    function getMesh2D()
    end
    function transformTri2Quad()
    end
  end
end