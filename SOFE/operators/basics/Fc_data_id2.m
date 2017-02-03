classdef Fc_data_id2 < Fc_Data_Id % (f, v)
  methods % constructor
    function obj = Fc_data_id2(coeff, feSpace, codim, varargin) % [loc]
      obj = obj@Fc_Data_Id(coeff, feSpace, codim, varargin{:});
    end
  end
end
