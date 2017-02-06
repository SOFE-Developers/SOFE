classdef SalomeMesh < Mesh
  methods % constructor
    function obj = SalomeMesh(file, varargin) % [dimP]
      [nodes, elem] = SalomeMesh.import(['./meshes/library/' file '.dat']);
      obj = obj@Mesh(nodes, elem, varargin{:});
    end
  end
  methods(Static = true)
    function [nodes, elem] = import(file) % tet
      fprintf('LOAD SALOMEMESH ...');
      fid = fopen(file, 'rt');
      s = textscan(fid, '%s', 'delimiter', '\n');
      %
      nN = dlmread(file,'',[0 0 0 0]);
      nE = dlmread(file,'',[0 1 0 1]);
      nL = countLines(file);
      nodes       = dlmread(file,'',[1 1 nN 3]);
      if norm(nodes(:,3))==0, nodes = nodes(:,[1 2]); end
      elem    = dlmread(file,'',[nL-nE-1 2 nL-2 9]);
      nV = size(elem, 2);
      N = find(elem(:,nV)>0,1,'first');
      elem = elem(N:end,:);
      if nV==8, elem = elem(:,[1 2 4 3 5 6 8 7]); end
      fclose(fid);
      fprintf('DONE \n');
    end
  end
end