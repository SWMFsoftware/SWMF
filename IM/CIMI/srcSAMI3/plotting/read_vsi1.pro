

; ion velocity data (hydrogen)

  close,1

  vsi1   = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'vsiu1.dat',/f77_unformatted

  print,' reading vsi1 '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  vsi1(*,*,0:nl-1,*) = tmpt
  vsi1(*,*,nl,*)     = vsi1(*,*,0,*)

  close,1

  delvar,tmpr,etmp

  end
