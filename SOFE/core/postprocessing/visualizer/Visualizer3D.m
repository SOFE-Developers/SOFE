classdef Visualizer3D < Visualizer
  methods % constructor
    function obj = Visualizer3D(feSpace)
      obj = obj@Visualizer(feSpace);
    end
  end
  methods
    function h = patch(obj, U, varargin)
      obj.test(U);
      try n = varargin{1}.n; catch, n = 1; end
      try deform = varargin{1}.deform; catch, deform = false; end
      try I = varargin{1}.loc; catch, I = ':'; end
      I = obj.feSpace.mesh.topology.isSurface(I);
      [Y,X] = meshgrid(linspace(0,1,n+1));
      if obj.feSpace.element.isSimplex()
        idx = X(:)+Y(:)<=1;
        X = X(idx); Y = Y(idx);
      end
      if obj.feSpace.element.isSimplex()
        col = reshapeTop(n+1:-1:1, 1:(n+2)*(n+1)/2);
        col1 = col; col1(n+1:n:n^2+n+1) = 0; col1(col1 == 0) = [];
        col2 = col; col2(:,1) = []; col2(col2 == 0) = [];
        col3 = col; col3(1,:) = []; col3(:,1) = []; col3(col3==0) = [];
        col4 = col; col4(n+1:n:n^2+n+1) = 0; col4(1,:) = 0; col4(col4==0) = [];
        col1 = reshape(col1,[],1);
        col2 = reshape(col2,[],1);
        col3 = reshape(col3,[],1);
        col4 = reshape(col4,[],1);
        faces = [col1 col1+1 col2; col3 col3-1 col4]; %nESubxnVx1
        offset = permute((0:sum(I)-1)*(n+1)*(n+2)/2, [1 3 2]); % 1x1xnE
      else
        col1 = (1:n*(n+1))'; col1((n+1)*(1:n)) = [];
        faces = [col1 col1+1 col1+n+1 col1+n+2]; %nESubxnVx1
        faces = faces(:, [1 2 4 3]);
        offset = permute((0:sum(I)-1)*(n+1)^2, [1 3 2]); % 1x1xnE
      end
      faces = permute(bsxfun(@plus, faces, offset), [1 3 2]); % nESubxnExnV
      if deform
        value = obj.feSpace.evalDoFVector(U,[X(:) Y(:)],[],0,I); % nExnPxnW
        vertices = obj.feSpace.mesh.evalReferenceMap([X(:) Y(:)],0,I) + value; % nExnPxnW
        value = sum(reshape(permute(value, [2 1 3]), [], 3).^2, 2).^0.5; % (nP*nE)
        vertices = reshape(permute(vertices, [2 1 3]), [], 3); % (nP*nE)xnW
        h = patch('faces', reshape(faces, [], size(faces,3)), 'vertices', vertices, ... 
              'facevertexcdata',value,'facecolor','interp', ...
              'edgecolor','interp');
        axis tight, axis equal
      else
        value = obj.feSpace.evalDoFVector(U,[X(:) Y(:)], [], 0, I); % nExnVxnC
        if size(value,3)>1
          value = sum(value.^2,3).^0.5; % nVxnE
        end
        value = reshape(value', [], 1); % (nP*nE)
        vertices = obj.feSpace.mesh.evalReferenceMap([X(:) Y(:)], 0, I); % nExnVx3
        vertices = reshape(permute(vertices, [2 1 3]), [], 3); % (nV*nE)x3
        h = patch('faces', reshape(faces, [], size(faces,3)), 'vertices', vertices, ... 
              'facevertexcdata',value,'facecolor','interp', ...
              'edgecolor','interp');
        axis tight, axis equal
      end
    end
    function h = surf(obj, U, varargin)
      try N = varargin{1}.N; catch, N = 200; end
      try deform = varargin{1}.deform; catch, deform = false; end
      try map = varargin{1}.map; catch, map = @(u,v)[u,v,0*u]; end
      try curl = varargin{1}.curl; catch, curl = 0; end
      try grad = varargin{1}.grad; catch, grad = 0; end
      try factor = varargin{1}.factor; catch, factor = []; end
      if numel(N) == 1, N = N*ones(2,1); end
      [A,B] = meshgrid(linspace(0,1,N(1)), linspace(0,1,N(2)));
      P = map(A(:),B(:)); % nPx3
      if curl 
        W = obj.feSpace.evalDoFVector(U,{P}, [], 1); % nPxnCxnD
        C = zeros(size(W(:,:,1)));
        C(:,1) = W(:,3,2) - W(:,2,3);
        C(:,2) = W(:,1,3) - W(:,3,1);
        C(:,3) = W(:,2,1) - W(:,1,2);
        W = C;
      elseif grad
        W = permute(obj.feSpace.evalDoFVector(U,{P}, [], 1), [1 3 2]); % nPxnCxnD
      else
        W = obj.feSpace.evalDoFVector(U,{P}, [], 0); % nPxnC
      end
      if ~isempty(factor), W = bsxfun(@times, W, factor(P)); end
      P = reshape(P, [size(A) 3]);
      if size(W,2) == 1
        W = reshape(W, size(A));
        h = surf(P(:,:,1), P(:,:,2), P(:,:,3), W);
      else
        W = reshape(W,size(A,1),size(A,2),[]); % nPaxnPbxnC
        if deform
          P = P+W;
          h = surf(P(:,:,1), P(:,:,2), P(:,:,3), sum(W.^2,3).^0.5); 
        else
          try scale = varargin{1}.scale; catch, scale = 1.0; end
          try width = varargin{1}.width; catch, width = 2; end
          try normalize = varargin{1}.normalize; catch, normalize = true; end
          try abs = varargin{1}.abs; catch, abs = true; end
          try vectors = varargin{1}.vectors; catch, vectors = true; end
          try n = varargin{1}.n; catch, n = 40; end
          if numel(n) == 1, n = n*ones(2,1); end
          absW = sum(W.^2, 3).^0.5; % nPaxnPb
          if abs
            h = surf(P(:,:,1),P(:,:,2),P(:,:,3),absW);
          end
          if vectors
            if normalize, W = bsxfun(@rdivide, W, absW); W(W==Inf) = 0; end
            filter = ceil(N(:)./n);
            P = P(1:filter(1):end, 1:filter(2):end,:);
            W = W(1:filter(1):end, 1:filter(2):end,:);
            hold on;
            h = quiver3(P(:,:,1),P(:,:,2),P(:,:,3),W(:,:,1),W(:,:,2),W(:,:,3), scale, ...
                        'linewidth', width, 'color',[0 0 1]);
            hold off;
          end
        end
      end
      axis equal;
      shading interp;
    end
    function h = scatter(obj, U, varargin)
      obj.test(U);
      isT = obj.feSpace.element.isSimplex();
      try I = varargin{1}.loc; catch, I = ':'; end
      I = obj.feSpace.mesh.topology.isSurface(I);
      try
        N = varargin{1}.N/sqrt(sum(I)*0.5^isT);
      catch
        N = max(2,ceil(300/sqrt(sum(I)*0.5^isT)));
      end
      try deform = varargin{1}.deform; catch, deform = false; end
      [Px, Py] = meshgrid(linspace(0,1,N)');
      P = [Px(:) Py(:)];
      if isT
        P = P(sum(P,2)<=1,:);
      end
      Z = obj.feSpace.evalDoFVector(U, P, [], 0, I); % nExnPxnC
      P = obj.feSpace.mesh.evalReferenceMap(P, 0, I); % nExnPxnW
      if size(Z,3) == 1
        P = reshape(P, [], size(P,3));
        h = plot3k([P(:,1),P(:,2),P(:,3)], 'ColorData', Z(:));
      else
        if deform
          P = reshape(P + Z, [], size(P,3)); % (nE*nP)xnW
          Z = sum(Z.^2, 3).^0.5; % nExnP
          h = plot3k([P(:,1),P(:,2), P(:,3)], 'ColorData', Z(:));
        else
          try scale = varargin{1}.scale; catch, scale = 1.0; end
          try width = varargin{1}.width; catch, width = 4; end
          try n = varargin{1}.n/sqrt(sum(I)*0.5^isT); catch, n = 50/sqrt(sum(I)*0.5^isT); end
          absZ = sum(Z.^2, 3).^0.5;
          Z = bsxfun(@rdivide, Z, absZ); Z(Z==Inf) = 0;
          h = plot3k(reshape(P,[],3), 'ColorData', absZ(:));
          P = reshape(P(:,1:ceil(N/n)^2:end,:), [], size(P,3)); % (nE*nP)xnW
          Z = reshape(Z(:,1:ceil(N/n)^2:end,:), [], size(Z,3)); % (nE*nP)xnW
          hold on
          quiver3(P(:,1),P(:,2),P(:,3),Z(:,1),Z(:,2),Z(:,3), scale, 'linewidth', width);
          hold off
        end
      end
      axis equal; axis tight;
    end
    function h = surfFH(obj, F, varargin)      
      try N = varargin{1}.N; catch, N = 200; end
      try map = varargin{1}.map; catch, map = @(u,v)[u,v,0*u]; end
      if numel(N) == 1, N = N*ones(2,1); end
      [A,B] = meshgrid(linspace(0,1,N(1)), linspace(0,1,N(2)));
      P = map(A(:),B(:)); % nPx3
      W = reshape(F(P), size(A)); % nPxnC
      P = reshape(P, [size(A) 3]);
      h = surf(P(:,:,1), P(:,:,2), P(:,:,3), W); shading interp;
      gs = obj.feSpace.mesh.topology.getGlobalSearcher();
      diam = gs.diam';
      axis(diam(:)); axis equal;
    end
  end
end