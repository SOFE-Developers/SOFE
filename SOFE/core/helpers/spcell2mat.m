function R = spcell2mat(C)
  sz(:,:,1) = cellfun(@(a)size(a,1),C);
  sz(:,:,2) = cellfun(@(a)size(a,2),C);
  assert(norm(max(sz(:,:,1),[],2) - min(sz(:,:,1),[],2))==0, 'Dimension mismatch');
  assert(norm(max(sz(:,:,2),[],1) - min(sz(:,:,2),[],1))==0, 'Dimension mismatch');
  %
  szI = [0; cumsum(sz(:,1,1))]; nI = numel(szI) - 1;
  szJ = [0; permute(cumsum(sz(1,:,2)), [2 1 3])]; nJ = numel(szJ) - 1;
  var = 2;
  switch var
    case 1
      R = sparse(szI(end), szJ(end));
      for k = 1:nI
        for l = 1:nJ
          R = R + ...
             [sparse(szI(end), szJ(l)), ...
             [sparse(szI(k), szJ(l+1)-szJ(l)); ...
              C{k,l}; ...
              sparse(szI(end)-szI(k+1), szJ(l+1)-szJ(l))], ...
              sparse(szI(end), szJ(end)-szJ(l+1))];
        end
      end
    case 2
      II = cell(nI,nJ);
      JJ = cell(nI,nJ);
      EE = cell(nI,nJ);
      for k = 1:nI
        for l = 1:nJ
          [II{k,l},JJ{k,l},EE{k,l}] = find(C{k,l});
          II{k,l} = II{k,l} + szI(k);
          JJ{k,l} = JJ{k,l} + szJ(l);
        end
      end
      R = sparse(cell2mat(II(:)), cell2mat(JJ(:)), cell2mat(EE(:)), szI(end), szJ(end));
  end
end