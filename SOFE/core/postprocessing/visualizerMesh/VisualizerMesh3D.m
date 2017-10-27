classdef VisualizerMesh3D < VisualizerMesh
  methods % constructor
    function obj = VisualizerMesh3D(mesh)
      obj = obj@VisualizerMesh(mesh);
    end
  end
  methods % display
    function show(obj, varargin) % [loc]
      c = caxis();
      fc = obj.mesh.topology.getEntity(2);
      Ib = obj.mesh.isBoundary(varargin{:});
      Is = obj.mesh.isSurface(varargin{:}) & ~Ib;
      h = trimesh(fc(Is,:), obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), obj.mesh.nodes(:,3));
      set(h,'facecolor',[0.5 0.7 0.2],'edgecolor','k');
      hold on
      if obj.mesh.topology.isSimplex, I = [1 2 3]; else, I = [1 2 4 3]; end
      h = trimesh(fc(Ib,I), obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), obj.mesh.nodes(:,3));
      hold off
      set(h,'facecolor',[0.5 0.8 0.5],'edgecolor','k');
      axis equal, axis tight, caxis(c);
    end
    function showInner(obj)
      nodes = obj.mesh.nodes;
      elem = obj.mesh.topology.getEntity(2);
      if obj.mesh.topology.isSimplex, I = [1 2 3]; else, I = [1 2 4 3]; end
      h = trimesh(elem(:,I), nodes(:,1), nodes(:,2), nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
    end
    function showCells(obj, I, flag)
      e2F = obj.mesh.getElem2Face();
      e2Ed = obj.mesh.getElem2Edge();
      elem = obj.mesh.topology.getEntity(3);
      II = cell(4,1);
      II{4} = I;
      II{3} = unique(e2F(I,:));
      II{2} = unique(e2Ed(I,:));
      II{1} = unique(elem(I,:));
      face = obj.mesh.topology.getEntity(2);
      nodes = obj.mesh.nodes;
      if obj.mesh.topology.isSimplex, I = [1 2 3]; else, I = [1 2 4 3]; end
      h = trimesh(face(II{3},I),nodes(:,1),nodes(:,2),nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
      if flag
        for i = 1:numel(flag)
          dim = str2double(flag(i));
          obj.showEntity(dim, II{dim+1});
        end
      end
    end
    function showEntity(obj, dim, varargin)
      if nargin < 3, I = (1:obj.mesh.getNumber(dim))'; else, I = varargin{1}; end
      center = obj.mesh.getCenter(dim);
      switch dim
        case 3
          color = [1 0 1];
        case 2
          color = [0 0 0];
        case 1
          color = [0 0.5 1];
        case 0
          color = [0 1 0.5];
      end
      text(center(I,1), center(I,2), center(I,3), num2str(I(:)),'Color',color,'FontSize', 18);
    end
  end
end