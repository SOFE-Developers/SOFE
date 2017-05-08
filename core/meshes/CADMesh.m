classdef CADMesh < Mesh
  properties
    type
  end
  methods % constructor
    function obj = CADMesh(file, varargin) % [dimP]
      obj = obj@Mesh([], []);
      if strcmp(file(end), 't')
        [nodes, elem] = CADMesh.importDATFormat(file);
      else
        [nodes, elem, obj.type] = CADMesh.importMESHFormat(file);
      end
      obj.initMesh(nodes, elem, varargin{:});
    end
  end
  methods(Static = true)
    function [nodes, elem] = importDATFormat(file)
      fprintf('LOAD CADMESH ...');
      fid = fopen(file, 'rt');
      %
      nN = dlmread(file,'',[0 0 0 0]);
      nE = dlmread(file,'',[0 1 0 1]);
      nL = countLines(file);
      nodes       = dlmread(file,'',[1 1 nN 3]);
      if norm(nodes(:,3))==0, nodes = nodes(:,[1 2]); end
      elem = dlmread(file,' ',nL-nE-1,2); elem = elem(:,1:end-1);
      nV = size(elem, 2);
      N = find(elem(:,nV)>0,1,'first');
      elem = elem(N:end,:);
      if nV==8, elem = elem(:,[1 2 4 3 5 6 8 7]); end
      fclose(fid);
      fprintf('DONE \n');
    end
    function [nodes, elem, type] = importMESHFormat(file)
      fprintf('LOAD CADMESH ...');
      fid = fopen(file, 'rt');
      if fid < 0, error('Cannot open mesh file'); end
      s = textscan(fid, '%s', 'delimiter', '\n');
      %
      iN = find(strcmp(s{1}, 'Vertices'), 1, 'first');
      nN = dlmread(file,'',[iN 0 iN 1]);
      nodes = dlmread(file,'',[iN+1 0 iN+nN 2]);
      %
      iT = find(strcmp(s{1}, 'Tetrahedra'), 1, 'first');
      nT = dlmread(file,'',[iT 0 iT 1]);
      elem = dlmread(file,'',[iT+1 0 iT+nT 4]);
      type = elem(:,5);
      elem = elem(:,1:4);
      fclose(fid);
      fprintf('DONE \n');
    end
  end
end