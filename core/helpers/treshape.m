function R = treshape(x, reverse)
  n = round( (-1 + sqrt(1+8*numel(x)))/2 );
  R = zeros(n);
  cnt = 1;
  if ~reverse
    for j = 1:n
      for i = 1:n
        if i+j>n+1, continue; end
        R(i,j) = cnt;
        cnt = cnt+1;
      end
    end
  else
    for i = 1:n
      for j = 1:n
        if i+j>n+1, continue; end
        R(i,j) = cnt;
        cnt = cnt+1;
      end
    end
  end
end
