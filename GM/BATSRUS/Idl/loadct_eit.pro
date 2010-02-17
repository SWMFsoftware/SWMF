pro loadct_eit, color
  ; load EIT table of given color: 'red', 'green', 'blue' or 'yellow'

  ; Find /Idl in the path
  ct_path = '.'
  list = strsplit(!path,':',/extract)
  for i = 0, n_elements(list)-1 do begin
     dir=list(i)
     if strpos(dir,'/Idl') gt strlen(dir)-6 then ct_path = dir+"/"
  endfor

  if n_elements(color) eq 0 then color=' '
  case color of
     'blue':  table=42
     'green': table=43
     'yellow':table=44
     'red':   table=45
     else:    table=41
  endcase

  loadct,file= ct_path+'color_table.eit',table
end


