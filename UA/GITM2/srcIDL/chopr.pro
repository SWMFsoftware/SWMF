function chopr, svalue, n
  if strlen(svalue) lt n then n = strlen(svalue)
  return, strmid(svalue, strlen(svalue)-n,n)
end
