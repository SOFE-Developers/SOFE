classdef SOFEClass < handle
  properties
    outputFlag = true;
  end
  methods % constructor
    function obj = SOFEClass()
    end
  end
  methods(Static = true)
    function R = getSOFEPath()
      R = pwd;
    end
    function R = getWorkPath()
      R = [pwd '/work'];
    end
    function R = getPluginPath()
      R = [pwd '/_plugins_'];
    end
    function install()
      addpath(genpath(SOFEClass.getPluginPath()));
      addpath(SOFEClass.getWorkPath());
      addpath(genpath([SOFEClass.getSOFEPath(),'/elements']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/feSpaces']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/helpers']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/meshes']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/operators']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/pdes']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/postprocessing']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/quadrature']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/solver']));
      more off
      openCase('demo');
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
end