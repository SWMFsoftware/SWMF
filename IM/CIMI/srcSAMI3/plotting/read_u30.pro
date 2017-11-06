

; interpolated u1p data

  u30   = fltarr(nx,ny,nl+1,nt)
  tmpt    = fltarr(nx,ny,nl,nt)
  tmp0    = fltarr(nx,ny,nl)

  close,1
  openr,1,dir+'u30.dat',/f77_unformatted
  print,' reading u30 '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor


  u30(*,*,0:nl-1,*) = tmpt
  u30(*,*,nl,*)     = u30(*,*,0,*)

  close,1

  delvar,tmpt,tmp0

  end
