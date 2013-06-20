pro loadct_extra, color
  ; load a color table specified by the 'color' argument string
  ; for '0' .. '40' use standard color table
  ; for other colors load color table from a file

  if n_elements(color) eq 0 then begin
     print, 'ERROR in loadct_extra: color argument is missing'
     return
  endif

  s = size(color)
  if s(0) ne 0 then begin
     print, 'ERROR in loadct_extra: color argument must be a scalar'
     return
  endif

  ; convert number value to a string if necessary
  if s(1) ne 7 then color = strtrim(string(fix(color)),2)

  if strlen(color) eq 0 then begin
     print, 'ERROR in loadct_extra: color has zero length'
     return
  endif

  ;  check if we match '0' ... '40'
  if stregex(color,'^[0-9]$',/boolean) or $
     stregex(color, '^[0-3][0-9]$', /boolean) or color eq '40' then begin
     loadct, fix(color)
     return
  endif

  ; Use this file for the extra color tables
  table1 = 'color_table_extra.ct'
  ; Find the color table in the path
  colortable1 = file_which(table1)
  if colortable1 eq '' then begin
     print, 'ERROR in loadct_extra: cannot find ',table1,' in the path=',!path
     return
  endif

  if n_elements(color) eq 0 then color=' '
  case color of
     'eityellow': loadct, file=colortable1, 41
     '171':   loadct, file=colortable1, 42
     '195':   loadct, file=colortable1, 43
     '284':   loadct, file=colortable1, 44
     '304':   loadct, file=colortable1, 45
      '94':   loadct, file=colortable1, 46
     '131':   loadct, file=colortable1, 47
     '193':   loadct, file=colortable1, 48
     '211':   loadct, file=colortable1, 49
     '335':   loadct, file=colortable1, 50
     '1600':  loadct, file=colortable1, 51
     '1700':  loadct, file=colortable1, 52
     '4500':  loadct, file=colortable1, 53
     else:    print, 'ERROR in loadct_extra: unknown table color=', color
  endcase

end


