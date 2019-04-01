function R = spLeftShiftNonZero(A)
  % shifts all non-zero entries to the left in each row
  assert(issparse(A), 'Matrix A must be sparse');
  [ii,~,ee] = find(A);
  u = unique([ii,ee],'rows');
  tmp = diff([0; find(diff(u(:,1))); size(u,1)]);
  R = reshapeTop(tmp, u(:,2));
%   R = mat2cell(u(:,2), tmp);
end



