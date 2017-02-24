function closeCase()
  name = isOpen();
  if ~isempty(name)
      folder = [SOFEClass.getWorkPath(),filesep,isOpen()];
      rmpath(genpath(folder));
  end
end