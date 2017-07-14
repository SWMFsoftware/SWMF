function tostr,value,nChars
  s = strcompress(string(long(value)),/remove_all)
  if (n_elements(nChars) eq 0) then nChars = strlen(s)
  if nChars lt strlen(s) then nChars = strlen(s)
  nZ = nChars-strlen(s)
  for i=0,nZ-1 do s = '0'+s
  return, s
end
