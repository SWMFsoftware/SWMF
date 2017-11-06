

; interpolated u1p data

  u1p0   = fltarr(nx,ny,nl+1,nt)
  tmpt    = fltarr(nx,ny,nl,nt)

  close,1
  openr,1,dir+'u1p0.dat',/f77_unformatted
  print,' reading u1p0 '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor


  u1p0(*,*,0:nl-1,*) = tmpt
  u1p0(*,*,nl,*)     = u1p0(*,*,0,*)

  close,1

  delvar,tmpt,tmp0

  end
