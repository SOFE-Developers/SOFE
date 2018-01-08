classdef Mass2 < PDE
  methods % constructor
    function obj = Mass2(C, FES)
      [m,n] = size(C);
      nnz = sum(cellfun(@(c)~isempty(c), FES));
      opList = cell(nnz,1);
      lhs.sys = cell(m,n); rhs.sys = cell(m, 1);
      k = 1;
      for i = 1:m
        for j = 1:n
          if isempty(C{i,j}), continue; end
          opList{k} = OpIdId(C{i,j}, 0, FES{j}, FES{i});
          lhs.sys{i,j} = {k};
          k = k+1;
        end
      end
      %
      obj = obj@PDE(opList, lhs, rhs);
    end
  end
end