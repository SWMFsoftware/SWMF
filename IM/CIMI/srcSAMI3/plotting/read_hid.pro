
; READ_HIC.PRO

; hidpv

  close,1
  hipcr    = fltarr(nf,nl,nt)
  hipcrtmp = fltarr(nf,nl)
  openr,1,dir+'hidpvu.dat',/f77_unformatted
  print,' reading hidpv'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,hipcrtmp
    hipcr(*,*,i) = hipcrtmp
  endfor

  hidpv             = fltarr(nf,nl+1,nt)
  hidpv(*,0:nl-1,*) = hipcr(*,*,*)
  hidpv(*,nl,*)     = hipcr(*,0,*)

  close,1

; hidphiv

  close,1
  hipcr    = fltarr(nf,nl,nt)
  hipcrtmp = fltarr(nf,nl)
  openr,1,dir+'hidphivu.dat',/f77_unformatted
  print,' reading hidphiv'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,hipcrtmp
    hipcr(*,*,i) = hipcrtmp
  endfor

  hidphiv             = fltarr(nf,nl+1,nt)
  hidphiv(*,0:nl-1,*) = hipcr(*,*,*)
  hidphiv(*,nl,*)     = hipcr(*,0,*)

  close,1

; hidpg

  close,1
  hipcr    = fltarr(nf,nl,nt)
  hipcrtmp = fltarr(nf,nl)
  openr,1,dir+'hidpgu.dat',/f77_unformatted
  print,' reading hidpg'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,hipcrtmp
    hipcr(*,*,i) = hipcrtmp
  endfor

  hidpg             = fltarr(nf,nl+1,nt)
  hidpg(*,0:nl-1,*) = hipcr(*,*,*)
  hidpg(*,nl,*)     = hipcr(*,0,*)

  close,1

; hidphig

  close,1
  hipcr    = fltarr(nf,nl,nt)
  hipcrtmp = fltarr(nf,nl)
  openr,1,dir+'hidphigu.dat',/f77_unformatted
  print,' reading hidphig'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,hipcrtmp
    hipcr(*,*,i) = hipcrtmp
  endfor

  hidphig             = fltarr(nf,nl+1,nt)
  hidphig(*,0:nl-1,*) = hipcr(*,*,*)
  hidphig(*,nl,*)     = hipcr(*,0,*)

  close,1


  end







