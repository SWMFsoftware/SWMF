
; READ_XYZ_S.PRO

; coordinate space

  close,1
  xr = fltarr(nz,nf,nl)
  openr,1,dir+'xsu.dat',/f77_unformatted
  print,' reading xs'
  readu,1,xr

  xs = fltarr(nz,nf,nl+1)
  xs(*,*,0:nl-1) = xr
  xs(*,*,nl)   = xs(*,*,0)

  close,1

  yr = fltarr(nz,nf,nl)
  openr,1,dir+'ysu.dat',/f77_unformatted
  print,' reading ys'
  readu,1,yr

  ys = fltarr(nz,nf,nl+1)
  ys(*,*,0:nl-1) = yr
  ys(*,*,nl)   = ys(*,*,0)

  close,1

  zr = fltarr(nz,nf,nl)
  openr,1,dir+'zsu.dat',/f77_unformatted
  print,' reading zs'
  readu,1,zr

  zs = fltarr(nz,nf,nl+1)
  zs(*,*,0:nl-1) = zr
  zs(*,*,nl)   = zs(*,*,0)

  close,1

  end
