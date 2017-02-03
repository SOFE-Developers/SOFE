function openCase(name)   
  closeCase();
  folder = [pwd,filesep,'work',filesep,name];
  addpath(genpath(folder));
  if isempty(isOpen())
    error('Work folder not found ');
  end
end