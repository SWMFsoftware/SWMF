PRO compute_axis, times, placement, btr, etr, curvar, xlab, xtl,	$
		basetime, moncheck, nl, ticktime, nminor

  ticktime = dblarr(20)
  xlab = strarr(20)
  timedum = intarr(6)

  nl = 1
  skip = 1

  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]

  pc = placement(curvar,0)

  if (times(pc,0,0) ge 65) then begin
    byear = 1900+times(pc,0,0)
    eyear = 1900+times(pc,1,0)
  endif else begin
    byear = 2000+times(pc,0,0)
    eyear = 2000+times(pc,1,0)
  endelse

  bmonth = times(pc,0,1)
  emonth = times(pc,1,1)
  bday = times(pc,0,2)
  eday = times(pc,1,2)
  bhour = times(pc,0,3)
  ehour = times(pc,1,3)
  bminute = times(pc,0,4)
  eminute = times(pc,1,4)
  bsecond = times(pc,0,5)
  esecond = times(pc,1,5)

  sbd = '0'+tostr(bday)
  sbd = strmid(sbd,strlen(sbd)-2,2)
  sed = '0'+tostr(eday)
  sed = strmid(sed,strlen(sed)-2,2)
  sbh = '0'+tostr(bhour)
  sbh = strmid(sbh,strlen(sbh)-2,2)
  seh = '0'+tostr(ehour)
  seh = strmid(seh,strlen(seh)-2,2)
  sbm = '0'+tostr(bminute)
  sbm = strmid(sbm,strlen(sbm)-2,2)
  sem = '0'+tostr(eminute)
  sem = strmid(sem,strlen(sem)-2,2)

  begt = times(pc,0,*)
  btr = double(0.0)
  c_a_to_r, begt, btr
  btr = btr - basetime(pc)

  endt = times(pc,1,*)
  etr = double(0.0)
  c_a_to_r, endt, etr
  etr = etr - basetime(pc)

  dt = etr - btr

  if (byear mod 4 eq 0) and (bmonth eq 2) then 		$
    dayofmon(1) = dayofmon(1)+1

  secperyear = 365.0*24.0*60.0*60.0
  secpermon  = float(dayofmon(bmonth-1))*24.0*60.0*60.0
  secperday  = 24.0*60.0*60.0
  secperhour = 60.0*60.0
  secpermin  = 60.0

  if dt ge 2.0*secperyear then begin

    step = (eyear - byear)/5 + 1
    step = 1
    nl = 0

    for j = byear, eyear+step, step do begin

      xlab(j-byear) = tostr(j)
      timedum = [j,1,0,0,0,0]
      c_a_to_r, timedum, dum
      ticktime(j-byear) = dum - basetime(pc)
      nl = nl + 1

    endfor

    xtl = 'UT Years'
    nminor = 6

  endif else begin

    if dt ge 2.0*secpermon then begin

      nmon = fix(dt/secpermon) - 1
      step = nmon/5 + 1
      nl = 0
      if byear eq eyear then exten = 0 else exten = 1

      for j=0, nmon+step, step do begin

	cmonth = bmonth + j
	cyear = byear
	if cmonth gt 12 then begin
	  if cmonth mod 12 ne 0 then begin
	    cyear = cyear + cmonth / 12
	    cmonth = cmonth mod 12
	  endif else begin
	    cyear = cyear + cmonth / 12 - 1
	    cmonth = 12
	  endelse
	endif

	xlab(nl) = strmid(moncheck,(cmonth-1)*3,3)
	if exten then xlab(nl) = xlab(nl) + ', ' + tostr(cyear)
	timedum = [cyear,cmonth,0,0,0,0]
	c_a_to_r, timedum, dum
	ticktime(nl) = dum - basetime(pc)
	nl = nl + 1

      endfor

      xlab(nl) = strmid(moncheck,(emonth-1)*3,3)
      nminor = 6
      if exten then xtl = 'UT Time'		$
      else xtl = tostr(byear)+' Universal Time'

    endif else begin

      if dt ge 2.0*secperday then begin

	nday = fix(dt/secperday) - 1
	step = nday/5 + 1
	nl = 0
	if bmonth eq emonth then exmon = 0 else exmon = 1
	if byear eq eyear then exyear = 0 else exyear = 1

	for j =	0, nday+step, step do begin

	  cday = bday + j
	  cmonth = bmonth
	  cyear = byear

; since we blindly incremented cday, we have to make sure that the right
; month and day are used. This is very easy to do, since the convering
; array to real subroutine does not care about too many days in a month - 
; so [91,2,29,0,0,0] = [91,3,1,0,0,0] ect.
; when we convert it back to an array, it is converted back normally.

	  timedum = [cyear,cmonth,cday,0,0,0]
	  c_a_to_r, timedum, dum
	  ticktime(nl) = dum - basetime(pc) 
	  c_r_to_a, timedum, dum

	  cmonth = timedum(1)
	  cday = timedum(2)

	  sday = '0'+tostr(cday)
	  sday = strmid(sday,strlen(sday)-2,2)
	  if exmon then 					$
	    xlab(nl) = strmid(moncheck,(cmonth-1)*3,3) + 	$
			' ' + sday				$
	  else xlab(nl) = sday

	  nl = nl + 1

	endfor

	nminor = step

	if exmon then begin

	  if exyear then					$
	    xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		  sbd + ', ' + tostr(byear) + ' to ' +	$
		  strmid(moncheck,(emonth-1)*3,3) + ' ' + 	$
		  sed + ', ' + tostr(eyear) + 		$
		  ' Universal Time'				$
	  else							$
	    xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		  sbd + ' to ' +			$
		  strmid(moncheck,(emonth-1)*3,3) + ' ' + 	$
		  sed + ', ' + tostr(eyear) + 		$
		  ' Universal Time'

	endif else						$
	  xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		sbd + ' to ' +				$
		sed + ', ' + tostr(eyear) + 		$
		' Universal Time'

      endif else begin

	if dt ge 2.0*secperhour then begin

	  nhour = fix(dt/secperhour) - 1
	  step = nhour/5 + 1
	  nl = 0
	  if bday eq eday then exday = 0 else exday = 1
	  if bmonth eq emonth then exmon = 0 else exmon = 1
	  if byear eq eyear then exyear = 0 else exyear = 1

	  for j = 0, nhour+step, step do begin

	    chour = bhour + j
	    cday = bday
	    cmonth = bmonth
	    cyear = byear

; since we blindly incremented cday, we have to make sure that the right
; month and day are used. This is very easy to do, since the convering
; array to real subroutine does not care about too many days in a month - 
; so [91,2,28,25,0,0] = [91,3,1,1,0,0] ect.
; when we convert it back to an array, it is converted back normally.

	    timedum = [cyear,cmonth,cday,chour,0,0]
	    c_a_to_r, timedum, dum
	    ticktime(nl) = dum - basetime(pc) 
	    c_r_to_a, timedum, dum

	    chour = timedum(3)

	    xlab(nl) = '0'+tostr(chour)
	    xlab(nl) = strmid(xlab(nl),strlen(xlab(nl))-2,2)

	    nl = nl + 1

	  endfor

	  nminor = step

	  if exday and (not exmon) and (not exyear) then 	$
	    if (bhour+bminute+bsecond eq 0) and			$
	       (ehour+eminute+esecond eq 0) and			$
               (bday eq eday-1) then exday = 0

	  if not exday then begin
	    xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		  sbd + ', ' + tostr(byear) + 		$
		  ' UT Hours'
	  endif else begin
	    if not exmon then begin
	      xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		    sbd + ' to ' +			$
		    sed + ', ' + tostr(eyear) + 	$
		    ' UT Hours'
	    endif else begin
	      if not exyear then begin
		xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		      sbd + ' to ' +			$
		      strmid(moncheck,(emonth-1)*3,3) + ' ' + 	$
		      sed + ', ' + tostr(eyear) + 	$
		      ' UT Hours'
	      endif else begin
		xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		      sbd + ', ' + tostr(byear) + 	$
		      ' to ' +					$
		      strmid(moncheck,(emonth-1)*3,3) + ' ' + 	$
		      sed + ', ' + tostr(eyear) + 	$
		      ' UT Hours'
	      endelse
	    endelse
	  endelse

	endif else begin

	  if dt ge 2.0*secpermin then do_min = 1 else do_min = 0

	  if do_min then ntotal = fix(dt/secpermin) - 1		$
	  else ntotal = fix(dt) - 1

	  step = ntotal/5 + 1

	  nl = 0
	  if bminute eq eminute then exmin = 0 else exmin = 1
	  if bhour eq ehour then exhour = 0 else exhour = 1
	  if bday eq eday then exday = 0 else exday = 1
	  if bmonth eq emonth then exmon = 0 else exmon = 1
	  if byear eq eyear then exyear = 0 else exyear = 1

	  for j = 0, ntotal+step, step do begin

	    if do_min then begin
	      csecond = 0
	      cminute = bminute + j
	    endif else begin
	      csecond = bsecond + j
	      cminute = bminute
	    endelse
	    chour = bhour
	    cday = bday
	    cmonth = bmonth
	    cyear = byear

; since we blindly incremented cday, we have to make sure that the right
; month and day are used. This is very easy to do, since the convering
; array to real subroutine does not care about too many days in a month - 
; so [91,2,28,25,0,0] = [91,3,1,1,0,0] ect.
; when we convert it back to an array, it is converted back normally.

	    timedum = [cyear,cmonth,cday,chour,cminute,csecond]
	    c_a_to_r, timedum, dum
	    ticktime(nl) = dum - basetime(pc) 
	    c_r_to_a, timedum, dum

	    if do_min then begin
              ctime = timedum(4) 
	      xlab(nl) = chopr('0'+tostr(timedum(3)),2)+chopr('0'+tostr(timedum(4)),2)+' UT'
            endif else begin
              ctime = timedum(5)
	    endelse

	    nl = nl + 1

	  endfor

	  nminor = step

	  if not exhour then begin
	      xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		    sbd + ', ' + tostr(byear) + 	$
		    ' ' + sbh + ':' + sbm
	  endif else begin
	    if not exday then begin
	        xtl = strmid(moncheck,(bmonth-1)*3,3) + ' ' + 	$
		      sbd + ', ' + tostr(byear) + 	$
		      ' ' + sbh + ':' + sbm +	$
		      ' to ' + seh + ':' + 		$
		      sem
	    endif else begin
	      if not exmon then begin
		  xtl = strmid(moncheck,(bmonth-1)*3,3) + 	$
			' ' + sbd + ' ' + 		$
			sbh + ':' + sbm +	$
			' to ' + 				$
			strmid(moncheck,(emonth-1)*3,3) + ' ' + $
			sed + ' ' + 		$
			seh + ':' + sem +	$
			', ' + tostr(eyear)
	      endif else begin
		if not exyear then begin
		    xtl = strmid(moncheck,(bmonth-1)*3,3) + 	$
			' ' + sbd + ' ' +		$
			sbh + ':' + sbm +	$
			' to ' + 				$
			strmid(moncheck,(emonth-1)*3,3) + ' ' + $
			sed + ' ' +			$
			seh + ':' + sem +	$
			', ' + tostr(eyear)
		endif else begin
		    xtl = strmid(moncheck,(bmonth-1)*3,3) + 	$
		      ' ' + sbd + ' ' +			$
		      sbh + ':' + sbm +	$
		      ', ' + tostr(byear) + 			$
		      ' to ' +					$
		      strmid(moncheck,(emonth-1)*3,3) + ' ' + 	$
		      sed + ' ' +			$
		      seh + ':' + sem +	$
		      ', ' + tostr(eyear)
		endelse
	      endelse
	    endelse
	  endelse
          if do_min then xtl = xtl + ' UT Minutes'		$
	  else xtl = xtl + ' UT Seconds'

	endelse

      endelse

    endelse

  endelse

  nl = nl - 1

;  loc = where(ticktime gt btr and ticktime lt etr, count)
;  if count gt 0 then begin
;    ticktime = ticktime(loc)
;    xlab = xlab(loc)
;    nl = count
;  endif

  if ticktime(0) lt btr then begin

    xlab = xlab(1:nl+1)
    ticktime = ticktime(1:nl)
    nl = nl - 1

  endif

  if (ticktime(nl) gt etr) or (ticktime(nl) lt btr) then begin

    xlab = [xlab(0:nl-1), xlab(nl+1)]
    ticktime = ticktime(0:nl-1)
    nl = nl - 1

  endif

  if nminor eq 1 then nminor = 4

  return

end
  
