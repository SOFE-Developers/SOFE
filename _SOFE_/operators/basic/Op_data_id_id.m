classdef Op_data_id_id < Op_data_Id_Id % ( c*u, v )
  methods % constructor
    function obj = Op_data_id_id(coeff, codim, fesTrial, varargin) % [fesTest loc]
      obj = obj@Op_data_Id_Id(coeff, codim, fesTrial, varargin{:});
    end
  end
end
