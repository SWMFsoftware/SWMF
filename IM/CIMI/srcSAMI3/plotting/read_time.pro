
; READ_TIME.PRO

  close,1

  time = fltarr(4,nt)
  openr,1,dir+'time.dat'
  print,' reading time'
  readf,1,time

  close,1

  print,'no of days?'
  read,days
  print,'days = ',fix(days)

  tm    = fltarr(nt)
  tm(*) = reform(time(1,*)) + reform(time(2,*))/60. + reform(time(3,*))/3600.
 
  for j = 1,days do begin
    print,'j = ',j
    for i = 1,nt-1 do begin
      if tm(i) lt tm(i-1) then tm(i) = tm(i) + 24.
     endfor
  endfor

  end
