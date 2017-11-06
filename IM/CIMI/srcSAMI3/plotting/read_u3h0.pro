

; interpolated u3h0 data

  u3h0   = fltarr(nx,ny,nl+1,nt)
  tmpt    = fltarr(nx,ny,nl,nt)
  tmp0    = fltarr(nx,ny,nl)

  close,1
  openr,1,dir+'u3h0.dat',/f77_unformatted
  print,' reading u3h0 '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  u3h0(*,*,0:nl-1,*) = tmpt
  u3h0(*,*,nl,*)     = u3h0(*,*,0,*)

  close,1

  delvar,tmpt,tmp0

  end
