

; interpolated electron density data and grid

;  nx   =  100
;  ny   =  100

  glat0   = fltarr(nx,ny,nl+1)
  glon0   = fltarr(nx,ny,nl+1)
  zalt0   = fltarr(nx,ny,nl+1)
  tmp0    = fltarr(nx,ny,nl)

  close,1
  openr,1,dir+'glat0.dat',/f77_unformatted
  print,' reading glat0 '
  readu,1,tmp0

  glat0(*,*,0:nl-1)  = tmp0
  glat0(*,*,nl)      = glat0(*,*,0)

  close,1
  openr,1,dir+'glon0.dat',/f77_unformatted
  print,' reading glon0 '
  readu,1,tmp0

  glon0(*,*,0:nl-1)  = tmp0
  glon0(*,*,nl)      = glon0(*,*,0)

  close,1
  openr,1,dir+'zalt0.dat',/f77_unformatted
  print,' reading zalt0 '
  readu,1,tmp0

  zalt0(*,*,0:nl-1)  = tmp0
  zalt0(*,*,nl)      = zalt0(*,*,0)

  vsi20   = fltarr(nx,ny,nl+1,nt)
  tmpt    = fltarr(nx,ny,nl,nt)

  close,1
  openr,1,dir+'vsi20.dat',/f77_unformatted
  print,' reading vsi20 '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor
;  readu,1,tmpt

  vsi20(*,*,0:nl-1,*) = tmpt
  vsi20(*,*,nl,*)     = vsi20(*,*,0,*)

  close,1

  delvar,tmpt,tmp0

  end
