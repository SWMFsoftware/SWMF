

; ion temperature data (hydrogen)

  close,1

  ti1   = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'tiu1.dat',/f77_unformatted

  print,' reading ti1 '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  ti1(*,*,0:nl-1,*) = tmpt
  ti1(*,*,nl,*)     = ti1(*,*,0,*)

  close,1

  delvar,tmpr,etmp

  end
