function R = padcell2mat(C)
  % concatenates (Nx1)-cellarray of [...]xnPx[...]x[...] along 1st dimension
  % by filling with NaN along 2nd dimension
  % 3rd & 4th dimension must match!
  N = numel(C);
  sizeVec = cellfun(@(x)size(x,2),C); 
  mx = max(sizeVec);
  R = cell(size(C));
  for k = 1:N
    R{k} = NaN(size(C{k},1), mx, size(C{k},3), size(C{k},4));
    R{k}(:,1:sizeVec(k),:,:) = C{k};
  end
  R = cell2mat(R);
end