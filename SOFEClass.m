classdef SOFEClass < handle
  properties(Constant)
    outputFlag = true;
    nBlockGS = 1;
    mem = 4; % [GBy]
  end
  methods % constructor
    function obj = SOFEClass()
    end
  end
  methods(Static = true) % path
    function R = getSOFEPath()
      R = [pwd '/_SOFE_'];
    end
    function R = getWorkPath()
      R = [pwd '/work'];
    end
    function R = getPluginPath()
      R = [pwd '/plugins'];
    end
    function unlock()
      addpath(genpath([SOFEClass.getSOFEPath(),'/elements']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/feSpaces']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/helpers']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/meshes']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/operators']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/pdes']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/postprocessing']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/quadrature']));
      addpath(genpath([SOFEClass.getSOFEPath(),'/solver']));
      SOFEClass.open('demo');
      more off
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
    function plugin(name)
      addpath(genpath([SOFEClass.getPluginPath(),'/', name]));
    end
    function plugout(varargin)
      if nargin > 0
        rmpath(genpath([SOFEClass.getPluginPath(),'/', varargin{:}]));
      else
        addpath(genpath(SOFEClass.getPluginPath()));
        rmpath(genpath(SOFEClass.getPluginPath()));
      end
    end
  end
  methods(Static = true)
    function R = getElementsPerBlock(nB, nQ, nC, nD)
      R = floor(1/3*10^9*SOFEClass.mem/8/nB^2/nQ/nC/nD);
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