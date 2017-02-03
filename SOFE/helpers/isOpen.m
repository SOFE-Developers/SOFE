function name = isOpen()
  x = path;
  i = strfind(x,'/work/');
  if ~isempty(i) 
      j = strfind(x(i(1):end),':');
      name = x(i(1)+6:i(1)+j(1)-2);
  else
      name = [];
  end
end