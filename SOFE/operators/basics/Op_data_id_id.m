classdef Op_data_id_id < Op_data_Id_Id % ( c*u, v )
  methods % constructor
    function obj = Op_data_id_id(coeff, codim, feSpaceTrial, varargin) % [feSpaceTest loc]
      obj = obj@Op_data_Id_Id(coeff, codim, feSpaceTrial, varargin{:});
    end
  end
end
