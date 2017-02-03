classdef Mass < PDE
  methods % constructor
    function obj = Mass(c, fes)
      if ~iscell(fes)
        fes = {fes};
      end
      nD = numel(fes);
      empty = cell(nD, 1);
      op = cell(nD, nD);
      if ~iscell(c), c = num2cell(c); end
      for i = 1:nD
        for j = 1:nD
          if ~isempty(c{i,j})
            op{i,j} = {Op_data_Id_Id(c{i,j}, 0, fes{i}, fes{j})};
          end
        end
      end
      obj = obj@PDE(op, empty);
    end
  end
end