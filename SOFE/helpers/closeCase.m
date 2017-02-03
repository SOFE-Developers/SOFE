function closeCase()
  name = isOpen();
  if ~isempty(name)
      folder = [pwd,filesep,'work',filesep,isOpen()];
      rmpath(genpath(folder));
  end
end