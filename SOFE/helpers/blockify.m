function R = blockify(matrixFamily)
  % matrixFamily ... NxNxI
  % R ... (I*N)x(I*N)
  [~,N,I] = size(matrixFamily);
  j = kron((1:I*N)', ones(N,1));
  i = permute(reshape(j,[N,N,I]), [2 1 3]);
  R = sparse(i(:),j(:),matrixFamily(:));
end