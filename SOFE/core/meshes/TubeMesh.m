classdef TubeMesh < RegularMesh
  methods % constructor
    function obj = TubeMesh(N, rad, len, varargin)
      N = [N ceil(len/(2*pi*rad)*N)];
      obj = obj@RegularMesh(N, [0 2*pi; 0 len], varargin{:});
      %
      map(:,1) = (1:N(1)+1:prod(N+1))';
      map(:,2) = map(:,1) + N(1);
      %
      iNVec = (1:obj.topology.getNumber(0))';
      iNVec(map(:,2)) = iNVec(map(:,1));
      [nodes, elem] = TorusMesh.removeNodes(obj.nodes, iNVec(obj.topology.getEntity(2)));
      obj.topology.update(elem);
      obj.nodes = [rad*cos(nodes(:,1)), rad*sin(nodes(:,1)), nodes(:,2)];
      obj.dimW = 3;
    end
  end
end