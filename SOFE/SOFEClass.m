classdef SOFEClass < handle
  properties
    outputFlag = true;
  end
  methods % constructor
    function obj = SOFEClass()
    end
  end
  methods
    function output(obj, str, type)
      if type <= obj.outputFlag
        breaker = '\n';
        switch type
          case 1
            fprintf(['#++++++++++++++++++++++++++', breaker]);
            fprintf(['#  ', str, '\n']);
            fprintf(['#++++++++++++++++++++++++++', breaker]);
          case 2
            fprintf(['# ', str, breaker]);
          case 3
            fprintf(['# --> ', str, breaker]);
          case 4
            fprintf([str, breaker]);
          otherwise
            fprintf(str);
        end
      end
    end
  end
  methods(Static = true)
    function R = getWorkPath()
      R = [pwd '/work'];
    end
    function R = getPluginPath()
      R = '';
    end
  end
end