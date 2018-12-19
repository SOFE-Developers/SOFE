classdef Visualizer2D < Visualizer
  methods % constructor
    function obj = Visualizer2D(feSpace)
      obj = obj@Visualizer(feSpace);
    end
  end
  methods
    function h = patch(obj, U, varargin)
      obj.test(U);
      try n = varargin{1}.n; catch, n = 1; end
      try deform = varargin{1}.deform; catch, deform = false; end
      delta = 1e-6; [Y,X] = meshgrid(linspace(delta,1-delta,n+1));
      if obj.feSpace.mesh.element.isSimplex()
        idx = X(:)+Y(:)<=1;
        X = X(idx); Y = Y(idx);
        Y = Y-eps; % circumvent duffy singular points
      end
      nE = obj.feSpace.mesh.topology.getNumber(2);
      if obj.feSpace.mesh.element.isSimplex()
        col = reshapeTop(n+1:-1:1, 1:(n+2)*(n+1)/2);
        col1 = col; col1(n+1:n:n^2+n+1) = 0; col1(col1 == 0) = [];
        col2 = col; col2(:,1) = []; col2(col2 == 0) = [];
        col3 = col; col3(1,:) = []; col3(:,1) = []; col3(col3==0) = [];
        col4 = col; col4(n+1:n:n^2+n+1) = 0; col4(1,:) = 0; col4(col4==0) = [];
        elem = [col1(:) col1(:)+1 col2(:); col3(:) col3(:)-1 col4(:)]; %nESubxnPx1
        offset = permute((0:nE-1)*(n+1)*(n+2)/2, [1 3 2]); % 1x1xnE
      else
        col1 = reshape(1:n*(n+1),n+1,n);
        col1(n+1,:) = []; col1 = col1(:);
        elem = [col1 col1+1 col1+n+1 col1+n+2]; %nESubxnPx1
        elem = elem(:,[1 2 4 3]);
        offset = permute((0:nE-1)*(n+1)^2, [1 3 2]); % 1x1xnE
      end
      elem = permute(bsxfun(@plus, elem, offset), [1 3 2]); % nESubxnExnP
      %
      if deform
        value = obj.feSpace.evalDoFVector(U,[X(:) Y(:)],[],0); % nExnPxnW
        vertices = obj.feSpace.mesh.evalReferenceMap([X(:) Y(:)],0) + value; % nExnPxnW
        value = sum(reshape(permute(value, [2 1 3]), [], 2).^2, 2).^0.5; % (nP*nE)
        vertices = reshape(permute(vertices, [2 1 3]), [], 2); % (nP*nE)xnW
        height = zeros(size(value));
        h = patch('faces', reshape(elem, [], size(elem,3)), 'vertices', [vertices height], ... 
                  'facevertexcdata',value,'facecolor','interp', ...
                  'edgecolor','interp');
        view(2), axis equal;
      else
        value = obj.feSpace.evalDoFVector(U,[X(:) Y(:)],[],0); % nExnPxnC
        if size(value,3)>1
          value = sum(value.^2,3).^0.5; % nExnP
        end
        value = reshape(value', [], 1); % (nP*nE)
        vertices = obj.feSpace.mesh.evalReferenceMap([X(:) Y(:)],0); % nExnPxnW
        vertices = reshape(permute(vertices, [2 1 3]), [], size(vertices,3)); % (nP*nE)x2
        if size(vertices,2)==2
%           vertices = [vertices, zeros(size(value))]; % height profile off
          vertices = [vertices, value]; % height profile on
        end
        h = patch('faces', reshape(elem, [], obj.feSpace.mesh.element.nV(end)), ...
                  'vertices', vertices , ... 
                  'facevertexcdata',value,'facecolor','interp', ...
                  'edgecolor','interp');
        view(2), axis equal; axis tight, grid off
      end
    end
    function h = surf(obj, U, varargin)
      obj.test(U);
      try N = varargin{1}.N; catch, N = 200; end
      if numel(N) == 1, N = N*ones(2,1); end
      try deform = varargin{1}.deform; catch, deform = false; end
      try curl = varargin{1}.curl; catch, curl = 0; end
      try factor = varargin{1}.factor; catch, factor = []; end
      try
        box = varargin{1}.box;
      catch
        gs = obj.feSpace.mesh.getGlobalSearcher();
        box = gs.diam';
      end
      [X,Y] = meshgrid(linspace(box(1), box(2), N(1)), ...
                       linspace(box(3), box(4), N(2)));
      P = [X(:) Y(:)];
      if curl 
        Z = obj.feSpace.evalDoFVector(U,{P}, [], 1); % nPxnCxnD
        Z = Z(:,2,1) - Z(:,1,2);
      else
        Z = obj.feSpace.evalDoFVector(U,{P}, [], 0); % nPxnC
        try, Z = Z(:,varargin{1}.component); catch, end
      end
      if ~isempty(factor), Z = bsxfun(@times, Z, factor(P)); end
      if size(Z,2) == 1
        h = surf(X,Y,reshape(Z,size(X))); shading interp
      else
        if deform
          absZ = reshape(sum(Z.^2, 2).^0.5, size(X));
          h = surf(X+reshape(Z(:,1),size(X)), Y+reshape(Z(:,2),size(X)), absZ); shading interp;
        else
          try scale = varargin{1}.scale; catch, scale = 1.0; end
          try width = varargin{1}.width; catch, width = 1; end
          try normalize = varargin{1}.normalize; catch, normalize = true; end
          try vectors = varargin{1}.vectors; catch, vectors = true; end
          try abs = varargin{1}.abs; catch, abs = true; end
          try n = varargin{1}.n; catch, n = 40; end
          if numel(n) == 1, n = n*ones(2,1); end
          absZ = sum(Z.^2, 2).^0.5;
          if abs
            surf(X,Y,0*X,reshape(absZ,size(X))); shading interp
          end
          if vectors
            if normalize, Z = bsxfun(@rdivide, Z, absZ); Z(Z==Inf) = 0; end
            try color = varargin{1}.color; catch, color = 'b'; end
            Z = reshape(Z, size(X,1), size(X,2), []);
            filter = ceil(N./n);
            X = X(1:filter(1):end, 1:filter(2):end);
            Y = Y(1:filter(1):end, 1:filter(2):end);
            Z = Z(1:filter(1):end, 1:filter(2):end,:);
            hold on
            h = quiver(X,Y,Z(:,:,1),Z(:,:,2), scale, 'linewidth', width,'color',color);
            hold off
          else
            h = surf(X,Y,reshape(absZ,size(X))); shading interp
          end
        end
      end
      view(2), axis equal; axis tight, grid off
    end
    function h = scatter(obj, U, varargin)
      obj.test(U);
      isT = obj.feSpace.element.isSimplex();
      nN = obj.feSpace.mesh.topology.getNumber(2);
      try
        N = max(2,varargin{1}.N/sqrt(nN*0.5^isT));
      catch
        N = max(2,ceil(300/sqrt(nN*0.5^isT)));
      end
      try codim = varargin{1}.codim; catch, codim = 0; end
      try deform = varargin{1}.deform; catch, deform = false; end
      points = linspace(0,1,N)';
      if codim == 0
        [pointsX, pointsY] = meshgrid(points);
        points = [pointsX(:) pointsY(:)];
      end
      if isT
        points = points(sum(points,2)<=1,:);
        points = points-eps; % circumvent duffy singular points
      end
      P = obj.feSpace.mesh.evalReferenceMap(points, 0); % nExnPxnW
      Z = obj.feSpace.evalDoFVector(U, points, [], 0); % nExnPxnC      
      if size(Z,3) == 1
        P = reshape(P, [], size(P,3));
        if size(P,2) == 2
          h = plot3k([P(:,1),P(:,2),Z(:)]);
        else
          h = plot3k([P(:,1),P(:,2), P(:,3)], 'ColorData', Z(:));
        end
      else
        if deform
          P = reshape(P + Z, [], size(P,3)); % (nE*nP)xnW
          Z = sum(Z.^2,3).^0.5;
          plot3k([P(:,1), P(:,2), Z(:)]);
        else
          try scale = varargin{1}.scale; catch, scale = 1.0; end
          try width = varargin{1}.width; catch, width = 1; end
          try normalize = varargin{1}.normalize; catch, normalize = true; end
          try vectors = varargin{1}.vectors; catch, vectors = 1; end
          try abs = varargin{1}.abs; catch, abs = true; end
          try n = varargin{1}.n; catch, n = 10; end
          absZ = sum(Z.^2, 3).^0.5;
          h = [];
          if abs
            if size(P,3) == 2
              h = plot3k([reshape(P,[],size(P,3)), zeros(size(P,1)*size(P,2),1)], 'ColorData', absZ(:));
            else
              h = plot3k(reshape(P,[],size(P,3)), 'ColorData', absZ(:));
            end
          end
          if vectors==1
            if normalize, Z = bsxfun(@rdivide, Z, absZ); Z(Z==Inf) = 0; end
            try color = varargin{1}.color; catch, color = 'b'; end
            P = reshape(P(:,1:n:end,:), [], size(P,3)); % (nE*nP)xnW
            Z = reshape(Z(:,1:n:end,:), [], size(Z,3)); % (nE*nP)xnW
            hold on
            if size(Z,2)==2
              quiver(P(:,1),P(:,2),Z(:,1),Z(:,2), scale, 'linewidth', width, 'color',color);
            else
              quiver3(P(:,1),P(:,2),P(:,3),Z(:,1),Z(:,2),Z(:,3), scale, 'linewidth', width);
            end
            hold off
          elseif vectors==2
            P = reshape(P(:,1:n:end,:), [], size(P,3)); % (nE*nP)xnW
            Z = reshape(Z(:,1:n:end,:), [], size(Z,3)); % (nE*nP)xnW
            hold on
            if size(Z,2)==4
              quiver(P(:,1),P(:,2),Z(:,1),Z(:,2), scale, 'linewidth', width, 'color','r');
              quiver(P(:,1),P(:,2),Z(:,3),Z(:,4), scale, 'linewidth', width, 'color','g');
            else
              quiver3(P(:,1),P(:,2),P(:,3),Z(:,1),Z(:,2),Z(:,3), scale, 'linewidth', width, 'color','r');
              quiver3(P(:,1),P(:,2),P(:,3),Z(:,4),Z(:,5),Z(:,6), scale, 'linewidth', width, 'color','g');
            end
            hold off
          end
        end
      end
      view(2), axis equal; axis tight, grid off
    end
    function h = surfFH(obj, F, varargin)
      try N = varargin{1}.N; catch, N = 200; end
      if numel(N) == 1, N = N*ones(2,1); end
      box = obj.feSpace.mesh.getDiam();
      [X,Y] = meshgrid(linspace(box(1), box(2), N(1)), ...
                       linspace(box(3), box(4), N(2)));
      Z = F([X(:) Y(:)]); % nPx1
      h = surf(X,Y,reshape(Z,size(X))); shading interp
      axis(box(:)); view(2), axis equal; axis tight
    end
  end
  methods(Static=true)
    function surfFunction(f, varargin) % [xgrid, ygrid]
      try
        xgrid = varargin{1};
      catch
        xgrid = linspace(0,1,400);
      end
      try
        ygrid = varargin{2};
      catch
        ygrid = linspace(0,1,400);
      end
      [X,Y] = meshgrid(xgrid, ygrid);
      surf(X,Y,reshape(f([X(:),Y(:)]),size(X)));
      shading interp
    end
  end
end