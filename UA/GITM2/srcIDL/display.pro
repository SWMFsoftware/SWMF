
pro display, vars

  nVars = n_elements(vars)

  if (nVars eq 0) then return

  nchop = floor(alog10(nVars))+1

  for iVar = 0,nVars-1 do $
    print, chopr('0000'+tostr(iVar),nchop),'. ',vars(iVar)

end
