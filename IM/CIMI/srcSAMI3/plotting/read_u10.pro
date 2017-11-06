

; interpolated u1p data

  u10   = fltarr(nx,ny,nl+1,nt)
  tmpt    = fltarr(nx,ny,nl,nt)

  close,1
  openr,1,dir+'u10.dat',/f77_unformatted
  print,' reading u10 '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor


  u10(*,*,0:nl-1,*) = tmpt
  u10(*,*,nl,*)     = u10(*,*,0,*)

  close,1

  delvar,tmpt,tmp0

  end
