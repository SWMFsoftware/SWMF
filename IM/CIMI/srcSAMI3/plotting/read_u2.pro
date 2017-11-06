

; u2

  close,1

  u2     = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'u2u.dat',/f77_unformatted

  print,' reading dene '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  u2(*,*,0:nl-1,*) = tmpt
  u2(*,*,nl,*)     = u2(*,*,0,*)

  close,1

  delvar,tmpr,etmp

  end
