
; READ_INDICES_DIR.PRO

; set the color table
  
  device,true=24
  device,de=0
  device,retain=2
  loadct,41

; sami3tgcm/primo

  nz    = 204
  nion  = 7
  nf    = 80
  nl    = 96

;  nz    = 101
;  nion  = 7
;  nf    = 124
;  nl    = 96

  print,'input nt:'
  read,nt

  nx = 200  ; latitude
  ny = 200  ; altitude

  nnx = nl+1
  nny = nf-1

dir='./'

; READ_TIME.PRO

  close,1

  time = dblarr(5,nt)
;  time = dblarr(4,nt)
  openr,1,dir+'time.dat'
  print,' reading time'
  readf,1,time

  close,1

  print,'no of days?'
  read,days
  print,'days = ',fix(days)

  tm    = dblarr(nt)
  tm(*) = reform(time(1,*)) + reform(time(2,*))/60. + reform(time(3,*))/3600.
 
  for j = 1,days do begin
    print,'j = ',j
    for i = 1,nt-1 do begin
      if tm(i) lt tm(i-1) then tm(i) = tm(i) + 24.
     endfor
  endfor


; READ_GEO_COORD.PRO

; magnetic longitude

  close,1
  glonr = dblarr(nz,nf,nl)
  openr,1,dir+'glonu.dat',/f77_unformatted
  print,' reading glon'
  readu,1,glonr

  glon = dblarr(nz,nf,nl+1)
  glon(*,*,0:nl-1) = glonr
  glon(*,*,nl)   = glon(*,*,0)

; magnetic latitude

  close,1
  glatr = dblarr(nz,nf,nl)
  openr,1,dir+'glatu.dat',/f77_unformatted
  print,' reading glat'
  readu,1,glatr

  glat = dblarr(nz,nf,nl+1)
  glat(*,*,0:nl-1) = glatr
  glat(*,*,nl)   = glat(*,*,0)

; altitude

  close,1
  zaltr = dblarr(nz,nf,nl)
  openr,1,dir+'zaltu.dat',/f77_unformatted
  print,' reading zalt'
  readu,1,zaltr

  zalt = dblarr(nz,nf,nl+1)
  zalt(*,*,0:nl-1) = zaltr
  zalt(*,*,nl)   = zalt(*,*,0)

; READ_MAG_COORD.PRO

; s magnetic longitude

  close,1
  blonr = dblarr(nz,nf,nl)
  openr,1,dir+'blonu.dat',/f77_unformatted
  print,' reading blon'
  readu,1,blonr

  blon = dblarr(nz,nf,nl+1)
  blon(*,*,0:nl-1) = blonr
  blon(*,*,nl)   = blon(*,*,0)

; s magnetic latitude

  close,1
  blatr = dblarr(nz,nf,nl)
  openr,1,dir+'blatu.dat',/f77_unformatted
  print,' reading blat'
  readu,1,blatr

  blat = dblarr(nz,nf,nl+1)
  blat(*,*,0:nl-1) = blatr
  blat(*,*,nl)   = blat(*,*,0)

; s altitude

  close,1
  baltr = dblarr(nz,nf,nl)
  openr,1,dir+'baltu.dat',/f77_unformatted
  print,' reading balt'
  readu,1,baltr

  balt = dblarr(nz,nf,nl+1)
  balt(*,*,0:nl-1) = baltr - 6370.
  balt(*,*,nl)   = balt(*,*,0)

; p magnetic longitude

  close,1
  blonr = dblarr(nz+1,nf+1,nl)
  openr,1,dir+'blonpu.dat',/f77_unformatted
  print,' reading blonp'
  readu,1,blonr

  blonp = dblarr(nz+1,nf+1,nl+1)
  blonp(*,*,0:nl-1) = blonr
  blonp(*,*,nl)   = blonp(*,*,0)

; p magnetic latitude

  close,1
  blatr = dblarr(nz+1,nf+1,nl)
  openr,1,dir+'blatpu.dat',/f77_unformatted
  print,' reading blatp'
  readu,1,blatr

  blatp = dblarr(nz+1,nf+1,nl+1)
  blatp(*,*,0:nl-1) = blatr
  blatp(*,*,nl)   = blatp(*,*,0)

; p altitude

  close,1
  baltr = dblarr(nz+1,nf+1,nl)
  openr,1,dir+'baltpu.dat',/f77_unformatted
  print,' reading baltp'
  readu,1,baltr

  baltp = dblarr(nz+1,nf+1,nl+1)
  baltp(*,*,0:nl-1) = baltr - 6370.
  baltp(*,*,nl)   = baltp(*,*,0)

; s coordinate space

  close,1
  xr = dblarr(nz,nf,nl)
  openr,1,dir+'xsu.dat',/f77_unformatted
  print,' reading xs'
  readu,1,xr

  xs = dblarr(nz,nf,nl+1)
  xs(*,*,0:nl-1) = xr
  xs(*,*,nl)   = xs(*,*,0)

  close,1

  yr = dblarr(nz,nf,nl)
  openr,1,dir+'ysu.dat',/f77_unformatted
  print,' reading ys'
  readu,1,yr

  ys = dblarr(nz,nf,nl+1)
  ys(*,*,0:nl-1) = yr
  ys(*,*,nl)   = ys(*,*,0)

  close,1

  zr = dblarr(nz,nf,nl)
  openr,1,dir+'zsu.dat',/f77_unformatted
  print,' reading zs'
  readu,1,zr

  zs = dblarr(nz,nf,nl+1)
  zs(*,*,0:nl-1) = zr
  zs(*,*,nl)   = zs(*,*,0)

  close,1

; p coordinate space

  close,1
  xr = dblarr(nz+1,nf+1,nl)
  openr,1,dir+'xpu.dat',/f77_unformatted
  print,' reading xp'
  readu,1,xr

  xp = dblarr(nz+1,nf+1,nl+1)
  xp(*,*,0:nl-1) = xr
  xp(*,*,nl)   = xp(*,*,0)

  close,1

  yr = dblarr(nz+1,nf+1,nl)
  openr,1,dir+'ypu.dat',/f77_unformatted
  print,' reading yp'
  readu,1,yr

  yp = dblarr(nz+1,nf+1,nl+1)
  yp(*,*,0:nl-1) = yr
  yp(*,*,nl)   = yp(*,*,0)

  close,1

  zr = dblarr(nz+1,nf+1,nl)
  openr,1,dir+'zpu.dat',/f77_unformatted
  print,' reading zp'
  readu,1,zr

  zp = dblarr(nz+1,nf+1,nl+1)
  zp(*,*,0:nl-1) = zr
  zp(*,*,nl)   = zp(*,*,0)

  close,1

  end







