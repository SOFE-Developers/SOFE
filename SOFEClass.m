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
      R = [pwd '/SOFE'];
    end
    function R = getWorkPath()
      R = [pwd '/SOFE/work'];
    end
    function R = getPluginPath()
      R = [pwd '/SOFE/_plugins_'];
    end
    function unlock()
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
      SOFEClass.open('demo');
    end
    function open(name)   
      SOFEClass.close();
      folder = [SOFEClass.getWorkPath(),filesep,name];
      addpath(genpath(folder));
      if isempty(SOFEClass.isOpen())
        error('Work folder not found ');
      end
    end
    function R = isOpen()
      x = path();
      i = strfind(x,'/work/');
      if ~isempty(i) 
          j = strfind(x(i(1):end),':');
          R = x(i(1)+6:i(1)+j(1)-2);
      else
          R = [];
      end
    end
    function close()
      name = SOFEClass.isOpen();
      if ~isempty(name)
          folder = [SOFEClass.getWorkPath(),filesep, SOFEClass.isOpen()];
          rmpath(genpath(folder));
      end
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