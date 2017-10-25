classdef VisualizerMeshTet < VisualizerMesh
  methods % constructor
    function obj = VisualizerMeshTet(mesh)
      obj = obj@VisualizerMesh(mesh);
    end
  end
  methods % display
    function show(obj, varargin) % [loc]
      top = obj.mesh.topology;
      c = caxis();
      fc = top.getEntity(2);
      Ib = obj.mesh.isBoundary(varargin{:});
      Is = obj.mesh.isSurface(varargin{:}) & ~Ib;
      h = trimesh(fc(Is,:), top.nodes(:,1), top.nodes(:,2), top.nodes(:,3));
      set(h,'facecolor',[0.5 0.7 0.2],'edgecolor','k');
      hold on
      h = trimesh(fc(Ib,:), top.nodes(:,1), top.nodes(:,2), top.nodes(:,3));
      hold off
      set(h,'facecolor',[0.5 0.8 0.5],'edgecolor','k');
      axis equal, axis tight, caxis(c);
    end
    function showInner(obj)
      top = obj.mesh.topology;
      nodes = top.getEntity(0);
      elem = top.getEntity(2);
      h = trimesh(elem, nodes(:,1), nodes(:,2), nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
    end
    function showCells(obj, I, flag)
      top = obj.mesh.topology;
      e2F = top.getElem2Face();
      e2Ed = top.getElem2Edge();
      elem = top.getEntity(3);
      II = cell(4,1);
      II{4} = I;
      II{3} = unique(e2F(I,:));
      II{2} = unique(e2Ed(I,:));
      II{1} = unique(elem(I,:));
      face = top.getEntity(2);
      nodes = top.getEntity(0);
      h = trimesh(face(II{3},:),nodes(:,1),nodes(:,2),nodes(:,3));
      set(h,'facecolor','none','edgecolor','k');
      if flag
        for i = 1:numel(flag)
          dim = str2double(flag(i));
          obj.showEntity(dim, II{dim+1});
        end
      end
    end
    function showEntity(obj, dim, varargin)
      top = obj.mesh.topology;
      if nargin < 3, I = (1:top.getNumber(dim))'; else, I = varargin{1}; end
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