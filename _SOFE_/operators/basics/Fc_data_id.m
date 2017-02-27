classdef Fc_data_id < Fc_Data_Id % (f, v)
  methods % constructor
    function obj = Fc_data_id(coeff, feSpace, codim, varargin) % [loc]
      obj = obj@Fc_Data_Id(coeff, feSpace, codim, varargin{:});
    end
  end
end
