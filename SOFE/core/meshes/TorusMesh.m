classdef TorusMesh < RegularMesh
  methods % constructor
    function obj = TorusMesh(N, R, r, isTri)
      obj = obj@RegularMesh(N, [0 2*pi; 0 2*pi], isTri);
      %
      maps = cell(obj.element.dimension,1);
      maps{1}(:,1) = (1:N(1)+1:prod(N+1))';
      maps{1}(:,2) = maps{1}(:,1) + N(1);
      maps{2}(:,1) = (1:N(1)+1)';
      maps{2}(:,2) = maps{2}(:,1) + prod(N) + N(2);
      %
      iNVec = (1:obj.topology.getNumber(0))';
      iNVec(maps{1}(:,2)) = iNVec(maps{1}(:,1));
      iNVec(maps{2}(:,2)) = iNVec(maps{2}(:,1));
      [nodes, elem] = TorusMesh.removeNodes(obj.nodes, iNVec(obj.topology.getEntity(2)));
      obj.topology.update(elem);
      obj.nodes = [(R+r*cos(nodes(:,1))).*cos(nodes(:,2)), ...
                   (R+r*cos(nodes(:,1))).*sin(nodes(:,2)), ...
                    r*sin(nodes(:,1))];
      obj.dimW = 3;
    end
  end
  methods(Static=true)
    function [nodes, elem] = removeNodes(nodes, elem)
      unode = unique(elem);
      nN = size(nodes,1); nNNew = numel(unode);
      nodes = nodes(unode,:);
      M = zeros(nN,1); M(unode) = (1:nNNew)';
      elem = M(elem);
    end
  end
end