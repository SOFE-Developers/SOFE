classdef Op_data_id_id2 < Op_data_Id_Id % ( c*u, v )
  methods % constructor
    function obj = Op_data_id_id2(coeff, codim, feSpaceTrial, varargin) % [feSpaceTest loc]
      obj = obj@Op_data_Id_Id(coeff, codim, feSpaceTrial, varargin{:});
    end
  end
end
