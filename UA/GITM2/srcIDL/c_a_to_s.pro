pro c_a_to_s, timearray, strtime, comp = comp, nice = nice, $
              IncludeDay=IncludeDay, IncludeSeconds=IncludeSeconds

  if (n_elements(comp) eq 0) then comp = 0
  if (n_elements(nice) eq 0) then nice = 0
  if (n_elements(IncludeDay) eq 0) then IncludeDay = 1
  if (n_elements(IncludeSeconds) eq 0) then IncludeSeconds = 1

  if (nice) then mon='JanFebMarAprMayJunJulAugSepOctNovDec' $
  else mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'

  sd = '0'+tostr(timearray(2))
  sd = strmid(sd,strlen(sd)-2,2)

  if (comp eq 0) then sm = strmid(mon,(timearray(1)-1)*3,3) $
  else begin
      sm = '0'+tostr(timearray(1))
      sm = strmid(sm,strlen(sm)-2,2)
  endelse

  year = timearray(0)
  if (nice) then begin
      if (year lt 65) then year = year + 2000
      if (year gt 65 and year lt 100) then $
        year = year + 1900
  endif else year = year mod 100
  sy = chopr('0'+tostr(year),2)
  sh = '0'+tostr(timearray(3))
  sh = strmid(sh,strlen(sh)-2,2)
  si = '0'+tostr(timearray(4))
  si = strmid(si,strlen(si)-2,2)
  ss = '0'+tostr(timearray(5))
  ss = strmid(ss,strlen(ss)-2,2)

  if (nice) then begin
      if (IncludeDay) then strtime = sm+'. '+sd+', '+sy+' ' else strtime=''
      strtime = strtime+sh+si
      if (IncludeSeconds) then strtime = strtime+':'+ss
      strtime = strtime+' UT'
  endif else begin
      if (comp eq 0) then $
        strtime = sd+'-'+sm+'-'+sy+' '+sh+':'+si+':'+ss+'.000' $
      else strtime = sy+sm+sd+sh+si+ss
  endelse

  RETURN

END

