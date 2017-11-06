
; READ_HIC.PRO

; hipcp

  close,1
  hipcr    = fltarr(nf,nl,nt)
  hipcrtmp = fltarr(nf,nl)
  openr,1,dir+'hipcpu.dat',/f77_unformatted
  print,' reading hipcp'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,hipcrtmp
    hipcr(*,*,i) = hipcrtmp
  endfor

  hipcp             = fltarr(nf,nl+1,nt)
  hipcp(*,0:nl-1,*) = hipcr(*,*,*)
  hipcp(*,nl,*)     = hipcr(*,0,*)

  close,1

; hipcphi

  close,1
  hipcr    = fltarr(nf,nl,nt)
  hipcrtmp = fltarr(nf,nl)
  openr,1,dir+'hipcphiu.dat',/f77_unformatted
  print,' reading hipcphi'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,hipcrtmp
    hipcr(*,*,i) = hipcrtmp
  endfor

  hipcphi             = fltarr(nf,nl+1,nt)
  hipcphi(*,0:nl-1,*) = hipcr(*,*,*)
  hipcphi(*,nl,*)     = hipcr(*,0,*)

  close,1

; hihcm

  close,1
  hipcr    = fltarr(nf,nl,nt)
  hipcrtmp = fltarr(nf,nl)
  openr,1,dir+'hihcmu.dat',/f77_unformatted
  print,' reading hihcm'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,hipcrtmp
    hipcr(*,*,i) = hipcrtmp
  endfor

  hihcm             = fltarr(nf,nl+1,nt)
  hihcm(*,0:nl-1,*) = hipcr(*,*,*)
  hihcm(*,nl,*)     = hipcr(*,0,*)

  close,1


  end







