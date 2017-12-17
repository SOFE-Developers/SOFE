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
      R = [pwd,filesep,'core'];
    end
    function R = getWorkPath()
      R = [pwd,filesep,'work'];
    end
    function R = getPluginPath()
      R = [pwd,filesep,'plugins'];
    end
    function unlock()
      addpath(genpath([SOFE.getCorePath(),filesep,'algorithm']));
			addpath(genpath([SOFE.getCorePath(),filesep,'elements']));
      addpath(genpath([SOFE.getCorePath(),filesep,'feSpaces']));
      addpath(genpath([SOFE.getCorePath(),filesep,'helpers']));
      addpath(genpath([SOFE.getCorePath(),filesep,'meshes']));
      addpath(genpath([SOFE.getCorePath(),filesep,'operators']));
      addpath(genpath([SOFE.getCorePath(),filesep,'pde']));
      addpath(genpath([SOFE.getCorePath(),filesep,'visualization']));
      addpath(genpath([SOFE.getCorePath(),filesep,'quadrature']));
      more off
    end
    function open(name)
      if isempty(strfind(path(),'_SOFE_')) %#ok<*STREMP>
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
      i = strfind(x,[filesep,'work',filesep]);
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
      addpath(genpath([SOFE.getPluginPath(),filesep, name]));
    end
    function plugout(varargin)
      if nargin > 0
        rmpath(genpath([SOFE.getPluginPath(),filesep, varargin{:}]));
      else
        addpath(genpath(SOFE.getPluginPath()));
        rmpath(genpath(SOFE.getPluginPath()));
      end
    end
    function pluggedin()
      X = path;
      I = strfind(X,[filesep,'plugins',filesep]);
      plugin = '';
      while ~isempty(I)
        x = X(I(1)+9:end);
        i = strfind(x,filesep);
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
      R = max(floor(1/3*10^9*SOFE.mem/8/nB^2/nQ/nC/nD), 1);
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
