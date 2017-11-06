

; interpolated u1p data

  u20   = fltarr(nx,ny,nl+1,nt)
  tmpt    = fltarr(nx,ny,nl,nt)
  tmp0    = fltarr(nx,ny,nl)

  close,1
  openr,1,dir+'u20.dat',/f77_unformatted
  print,' reading u20 '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor


  u20(*,*,0:nl-1,*) = tmpt
  u20(*,*,nl,*)     = u20(*,*,0,*)

  close,1

  delvar,tmpt,tmp0

  end
