function R = sortCol(M)
  % determine ascending sorting indices of each column
  % M must be integer valued
  [n,N] = size(M);
  J = kron((1:N)', ones(n,1));
  A = sparse(M(:), J,1)>0;
  out = reshapeTop(full(sum(A)));
  [ii,jj,~] = find(A);
  A = sparse(ii,jj,out(out>0));
  I = M + size(A,1)*(kron(1:N,ones(n,1))-1);
  R = full(A(I));
end