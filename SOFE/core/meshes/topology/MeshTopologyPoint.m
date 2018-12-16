classdef MeshTopologyPoint < MeshTopology
  methods % constructor
    function obj = MeshTopologyPoint(elem)
      obj = obj@MeshTopology(0);
      obj.update(elem);
      obj.nESub = 1;
      obj.nO = 0;
    end
    function update(obj, elem)
      obj.connectivity = cell(obj.dimP+1);
      obj.connectivity{obj.dimP+1,1} = elem;
      obj.connectivity{1,1} = (1:max(obj.connectivity{1,1}(:)))';
    end
  end
  methods % connectivity information
    function R = getOrientation(obj, dim, d, varargin) %#ok<INUSD> % [I]      
      R = [];
    end
  end  
end