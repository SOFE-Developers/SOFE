classdef Mass2 < PDE2
  methods % constructor
    function obj = Mass2(fes, varargin) % [N]
      if ~isempty(varargin), N = varargin{1}; else, N = 1; end
      opList = {OpIdId(1.0, 0, fes)};
      %      
      rhs.sys = cell(N, 1);
      lhs.sys = cell(N, N);
      for i = 1:N
        lhs.sys{i,i} = {1};
      end
      obj = obj@PDE2(opList, lhs, rhs);
    end
  end
end