classdef SOFE < handle
  properties(Constant)
    outputFlag = true;
    mem = 4; % [GBy]
  end
  methods % constructor
    function obj = SOFE()
    end
  end
  methods(Static = true) % path
    function R = getCorePath()
      R = [pwd '/core'];
    end
    function R = getWorkPath()
      R = [pwd '/work'];
    end
    function R = getPluginPath()
      R = [pwd '/plugins'];
    end
    function unlock()
      addpath(genpath([SOFE.getCorePath(),'/elements']));
      addpath(genpath([SOFE.getCorePath(),'/feSpaces']));
      addpath(genpath([SOFE.getCorePath(),'/helpers']));
      addpath(genpath([SOFE.getCorePath(),'/meshes']));
      addpath(genpath([SOFE.getCorePath(),'/operators']));
      addpath(genpath([SOFE.getCorePath(),'/pdes']));
      addpath(genpath([SOFE.getCorePath(),'/postprocessing']));
      addpath(genpath([SOFE.getCorePath(),'/quadrature']));
      addpath(genpath([SOFE.getCorePath(),'/solver']));
      more off
    end
    function open(name)
      if isempty(strfind(path(),'_SOFE_'))
        SOFE.unlock();
      end  
      SOFE.close();
      folder = [SOFE.getWorkPath(),filesep,name];
      addpath(folder);
      if isempty(SOFE.isOpen())
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
      name = SOFE.isOpen();
      if ~isempty(name)
          folder = [SOFE.getWorkPath(),filesep, SOFE.isOpen()];
          rmpath(genpath(folder));
      end
    end
    function plugin(name)
      addpath(genpath([SOFE.getPluginPath(),'/', name]));
    end
    function plugout(varargin)
      if nargin > 0
        rmpath(genpath([SOFE.getPluginPath(),'/', varargin{:}]));
      else
        addpath(genpath(SOFE.getPluginPath()));
        rmpath(genpath(SOFE.getPluginPath()));
      end
    end
    function pluggedin()
      X = path;
      I = strfind(X,'/plugins/');
      plugin = '';
      while ~isempty(I)
        x = X(I(1)+9:end);
        i = strfind(x,'/');
        if ~strcmp(plugin, x(1:i(1)-1))
          plugin = x(1:i(1)-2);
          fprintf([plugin '\n']);
        end
        I(1) = [];
      end
    end
  end
  methods(Static = true)
    function R = getElementsPerBlock(nB, nQ, nC, nD)
      R = floor(1/3*10^9*SOFE.mem/8/nB^2/nQ/nC/nD);
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