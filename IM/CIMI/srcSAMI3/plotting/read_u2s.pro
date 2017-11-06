

; E x B drift (bfield - u2s)

  close,1

  u2s  = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'u2su.dat',/f77_unformatted

  print,' reading u2s '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  u2s(*,*,0:nl-1,*) = tmpt
  u2s(*,*,nl,*)     = u2s(*,*,0,*)

  close,1
 
  u2smax = fltarr(nt)
  for i = 0,nt-1 do u2smax(i)=max(u2s(*,*,*,i))

  u2smin = fltarr(nt)
  for i = 0,nt-1 do u2smin(i)=min(u2s(*,*,*,i))

  delvar,tmpr,etmp

  end
