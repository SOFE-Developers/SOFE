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
      % Adaptively refines all cells in given subdomain
      %
      % Call:
      % m.adaptiveRefine(loc[,N])
      %
      % loc ... function handle, location of cells to be refined
      % N   ... 1x1, number of refinements, default N = 1
      %
    end
    function coarsen()
      % Coarsens all cells in given subdomain that have been refined
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
      % Rotates mesh by given angle around origin (2D only)
      %
      % Call:
      % m.rotate(alpha)
      %
      % alpha ... 1x1, angle in rad
      %
    end
    function scale()
      % Scales mesh by given vector
      %
      % Call:
      % m.scale(a)
      %
      % a ... 1x1 or nWx1, acaling vector
      %
    end
    function translate()
      % Shifts mesh by given translation vector
      %
      % Call:
      % m.translate(vec)
      % vec ... nWx1, translation vector
      %
    end
    function applyLinearMap()
      % Transforms mesh by linear map A
      %
      % Call:
      % m.applyLinearMap(A)
      %
      % A ... nWxnW, linear map
      %
    end
    % mesh information
    function getMeasure()
      % Computes measure of (selected) mesh entities of given dimension
      %
      % Call:
      % R = m.getMeasure(dim[, I])
      %
      % dim ... 1x1, dimension of entities
      % I   ... logical or integer index vector
      % R   ... nEx1, list of measures
      %
    end
    function getCenter()
      % Computes centroid of (selected) mesh entities of given dimension
      %
      % Call:
      % R = m.getCenter(dim[, I])
      %
      % dim ... 1x1, dimension of entities
      % I   ... logical or integer index vector
      % R   ... nExnW, list of coordinates 
      %
    end
    function getOuterNormal()
      % Computes outer unit  normal vector for boundary faces in local points
      % (2D only)
      %
      % Call:
      % R = m.getOuterNormal(points)
      %
      % points ... nPx1, local points
      % R   ... nExnPxnW, list of normal vectors in local points
      %
    end
    function getOuterNormal2()
      % Computes outer unit normal vector for boundary faces in langrange points
      % to given reference element
      %
      % Call:
      % R = m.getOuterNormal2(element)
      %
      % element ... Element, reference element
      % R   ... nDoFxnW, list of normal vectors in local points
      %                  (NaN if not on boundary)
      %
    end
    function getDiam()
      % Returns coordinates of bounding box to mesh
      %
      % Call:
      % R = m.getDiam()
      %
      % R   ... 2xnW, max and min for all coordinate directions
      %
    end
    function findEntity()
      % Returns if entity lies with at least one vertex in given subdomain
      %
      % Call:
      % R = m.findEntity(dim[, loc])
      %
      % dim ... 1x1, dimension of entity
      % loc ..., function handle, logical discription of subdomain
      % R   ... nEx1, logical vector
      %
    end
    function findEntityC()
      % Returns if entity lies with centroid in given subdomain
      %
      % Call:
      % R = m.findEntityC(dim[, loc])
      %
      % dim ... 1x1, dimension of entity
      % loc ..., function handle, logical discription of subdomain
      % R   ... nEx1, logical vector
      %
    end
    function isBoundary()
      % Returns if facet lies on boundary (AND in given subdomain)
      %
      % Call:
      % R = m.isBoundary([loc])
      %
      % loc ..., function handle, logical discription of subdomain
      % R   ... nEx1, logical vector
      %
    end
    function getBoundary()
      % Returns list of facets that lie on boundary (AND in given subdomain)
      %
      % Call:
      % R = m.getBoundary([loc])
      %
      % loc ..., function handle, logical discription of subdomain
      % R   ... nEx2, list of facets
      %
    end
    function isSurface()
      % Returns if facet lies on boundary (OF given subdomain)
      %
      % Call:
      % R = m.isSurface([loc])
      %
      % loc ..., function handle, logical discription of subdomain
      % R   ... nEx1, logical vector
      %
    end
    function getSurface()
      % Returns list of facets that lie on boundary (OF given subdomain)
      %
      % Call:
      % R = m.getSurface([loc])
      %
      % loc ..., function handle, logical discription of subdomain
      % R   ... nEx2, list of facets
      %
    end
    function isBoundaryNode()
      % Returns if node lies on boundary (AND given subdomain)
      %
      % Call:
      % R = m.isBoundaryNode([loc])
      %
      % loc ..., function handle, logical discription of subdomain
      % R   ... nEx1, logical vector
      %
    end
    function getBoundaryNode()
      % Returns coordinates of nodes that lie on boundary (AND given subdomain)
      %
      % Call:
      % R = m.getBoundaryNode([loc])
      %
      % loc ..., function handle, logical discription of subdomain
      % R   ... nExnW, list of coordinates
      %
    end
    % mesh visualization
    function show()
      % Visualizes mesh by means of the corresponding VisualizerMesh class
      %
      % Call:
      % m.show()
      %
    end
  end
end