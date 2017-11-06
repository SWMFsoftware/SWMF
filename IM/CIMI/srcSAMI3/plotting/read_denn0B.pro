

; interpolated neutral density data and grid

  nx   =  100
  ny   =  100

  blat0   = fltarr(nx,ny,nl+1)
  blon0   = fltarr(nx,ny,nl+1)
  balt0   = fltarr(nx,ny,nl+1)
  tmp0    = fltarr(nx,ny,nl)

  close,1
  openr,1,dir+'blat0.dat',/f77_unformatted
  print,' reading blat0 '
  readu,1,tmp0

  blat0(*,*,0:nl-1)  = tmp0
  blat0(*,*,nl)      = blat0(*,*,0)

  close,1
  openr,1,dir+'blon0.dat',/f77_unformatted
  print,' reading blon0 '
  readu,1,tmp0

  blon0(*,*,0:nl-1)  = tmp0
  blon0(*,*,nl)      = blon0(*,*,0)

  close,1
  openr,1,dir+'balt0.dat',/f77_unformatted
  print,' reading balt0 '
  readu,1,tmp0

  balt0(*,*,0:nl-1)  = tmp0
  balt0(*,*,nl)      = balt0(*,*,0)

  denn10B  = fltarr(nx,ny,nl+1,nt)
;  deni10B  = fltarr(nx,ny,nl+1,nt)
  tmpt    = fltarr(nx,ny,nl,nt)

  close,1
  openr,1,dir+'denn10B.dat',/f77_unformatted
  print,' reading denn10B '
;  openr,1,dir+'deni10B.dat',/f77_unformatted
;  print,' reading deni10B '
  readu,1,tmpt

  denn10B(*,*,0:nl-1,*) = tmpt
  denn10B(*,*,nl,*)     = denn10B(*,*,0,*)
;  deni10B(*,*,0:nl-1,*) = tmpt
;  deni10B(*,*,nl,*)     = deni10B(*,*,0,*)

  close,1

  delvar,tmpt,tmp0

  end
