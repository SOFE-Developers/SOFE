function R = leftShiftNonZero(A);
  % shifts all non-zero entries to the left in each row
  A = A.';
  R = zeros(size(A));
  I = (A == 0);
  R(~sort(I)) = A(~I);
  R = R';
end