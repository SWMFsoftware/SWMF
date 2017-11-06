

; electron temperature data

  close,1

  te   = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'teu.dat',/f77_unformatted

  print,' reading te '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  te(*,*,0:nl-1,*) = tmpt
  te(*,*,nl,*)     = te(*,*,0,*)

  close,1

  delvar,tmpr,etmp

  end
