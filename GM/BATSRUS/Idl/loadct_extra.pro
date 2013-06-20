pro loadct_extra, color
  ; load a color table specified by the 'color' argument string
  ; for '0' .. '40' use standard color table
  ; for other colors load color table from a file

  if n_elements(color) eq 0 then begin
     print, 'ERROR in loadct_extra: color argument is missing'
     return
  endif

  if strlen(color) eq 0 then begin
     print, 'ERROR in loadct_extra: color has zero length'
     return
  endif

  if strlen(color) eq 1 and $
     strlen(stregex(color,'[0-9]',/extract)) eq 1 then begin
     loadct, fix(color)
     return
  endif

  if strlen(color) eq 2 and $
     strlen(stregex(color, '[0-4][0-9]', /extract)) eq 2 then begin
     number = fix(color)
     if number le 40 then begin
        loadct, number
        return
     endif
  endif

  table1 = 'color_table_extra.ct'
  ; Find the color table in the path
  colortable1 = file_which(table1)
  if colortable1 eq '' then begin
     print, 'ERROR in loadct_extra: cannot find ',table1,' in the path=',!path
     return
  endif

  if n_elements(color) eq 0 then color=' '
  case color of
     '0':    loadct,0
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


