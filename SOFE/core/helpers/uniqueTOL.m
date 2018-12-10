function [z,ii,jj] = uniqueTOL(x, tol, varargin)
  % almost-uniqueness for SOFE by SFranz
  % [z,ii,jj] = uniquetol2(x,tol,varargin)
  %
  % tries to emulate Matlab's uniquetol in finding unique elements in x up to
  % specified tolerance tol
  % 
  % inspirated by uniquetol by Siyi Deng, which is flawed.
  % used characterisation of a point in n-dimensional space by n+1 distances
  % to different non-linear/planar/hyperplanar points and ordering in 1d
  %
  % SFranz, 2017
  %
  % PS: works for not to small tolerances (tol>>eps) and not to large points
  % ||x||<<1/eps, otherwise comparison with globar set of distance-points ist
  % numerical unstable. Alternative: use local distance points in each
  % sphere. This is not programmed, yet.
  %
  if size(x,1) == 1, x = x(:); end
  dim = size(x,2);
  if nargin < 2 || isempty(tol) || tol == 0
      [z,ii,jj] = unique(x,varargin{:});
      return;
  end
  [y,ii,jj] = unique(x,varargin{:});
  if dim==1
    isTol = [true;diff(y) > tol];
    z = y(isTol,:);
    bin = cumsum(isTol); % [n,bin] = histc(y,z);
    jj = bin(jj);
    ii = ii(isTol);
  else    
    N = size(y,1);                 % number of unique points
    d = zeros(N,dim+1);            % distance-field
    d(:,1) = sqrt(sum(y.^2,2));    % distance to zero
    for i=2:dim+1
      e = zeros(1,dim);e(i-1) = 1; % unit-vector
      d(:,i) = sqrt(sum((y-repmat(e,N,1)).^2,2)); % distance to unit-vector
    end
    for i=1:dim+1
      [~,ord] = sort(d(:,i));              % sort distances
      isTol = [true;diff(d(ord,i)) > tol]; % is distance-difference larger than tol?
                                           % not the distance to the first 
                                           % point in each cluster, but to 
                                           % the previous point
      d(ord,i) = cumsum(isTol);            % replace distance by sphere-number
    end
    [~,i,j] = unique(d,varargin{:});
    z = y(i,:);
    ii = ii(i);
    jj = j(jj);
  end
end