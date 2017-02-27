function openCase(name)   
  closeCase();
  folder = [SOFEClass.getWorkPath(),filesep,name];
  addpath(genpath(folder));
  if isempty(isOpen())
    error('Work folder not found ');
  end
end