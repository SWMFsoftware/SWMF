; READ-BASIC.PRO

; set the color table 
  
  device,true=24
  device,de=0
  device,retain=2
  loadct,41
  re = 6370.0

; sami3tgcm

;  nz    = 201
;  nion  = 7
;  nf    = 120
;  nl    = 96
;  nt    = 100

; kate zawdie

;  nz    = 201
;  nion  = 7
;  nf    = 200
;  nl    = 96
;  nt    = 14

; esf

;  nz = 101
;  nf = 202
;  nion = 7
;  nl = 96
;  nt = 101

; volland-stern

;  nz  = 160
;  nf  = 200
;  nl  =  90
;  nt  = 187

; weimer

  nz  = 204
  nf  = 124
  nl  =  96
  nt  =  1

  nx = 100  ; latitude
  ny = 200  ; altitude
  nr = 60   ; radius (for plasmasphere diagnostic)
  icorot = 1 ; set to 1 if potential includes coroation, otherwise 0

  nnx = nl+1
  nny = nf-1

;  dir0d1 = '/sdd1/pogo2/col/sami3_wvs/sami3_mpi-1.90cr_euvac93/2001/'
;  dir0h  = '/home/pogo2/col/sami3_wvs/sami3_mpi-1.90cr_euvac93vs/2001/'
;  dir= dir0b1 + 'day328/prerun/'
;  dir= dir0h + 'day33/run1_i2a/'

;  dir='/home/pogo2/col/sami3_wvs/sami3_mpi-1.90cr_euvac93vs/2001/day33/prerun_i2a/'

dir='/home/jack/huba/ionosphere/sami3/sami3_mpi-1.96/standard/code_newgrid/'

; READ_TIME.PRO

  close,1

;  time = fltarr(5,nt)
;  openr,1,dir+'time.dat'
;  print,' reading time'
;  readf,1,time

  close,1

; this version handles up to a week
;  print,'no of days?'
;  read,days
;  print,'days = ',fix(days)
  days = 7

  tm    = fltarr(nt)
  tm(*) = reform(time(1,*)) + reform(time(2,*))/60. + reform(time(3,*))/3600.
 
  for j = 1,days do begin
;    print,'j = ',j
    for i = 1,nt-1 do begin
      if tm(i) lt tm(i-1) then tm(i) = tm(i) + 24.
     endfor
  endfor


; READ_GEO_COORD.PRO

; magnetic longitude

  close,1
  glonr = fltarr(nz,nf,nl)
  openr,1,dir+'glonu.dat',/f77_unformatted
  print,' reading glon'
  readu,1,glonr

  glon = fltarr(nz,nf,nl+1)
  glon(*,*,0:nl-1) = glonr
  glon(*,*,nl)   = glon(*,*,0)

; magnetic(?) latitude

  close,1
  glatr = fltarr(nz,nf,nl)
  openr,1,dir+'glatu.dat',/f77_unformatted
  print,' reading glat'
  readu,1,glatr

  glat = fltarr(nz,nf,nl+1)
  glat(*,*,0:nl-1) = glatr
  glat(*,*,nl)   = glat(*,*,0)

; altitude

  close,1
  zaltr = fltarr(nz,nf,nl)
  openr,1,dir+'zaltu.dat',/f77_unformatted
  print,' reading zalt'
  readu,1,zaltr

  zalt = fltarr(nz,nf,nl+1)
  zalt(*,*,0:nl-1) = zaltr
  zalt(*,*,nl)   = zalt(*,*,0)

; READ_MAG_COORD.PRO

; magnetic longitude

  close,1
  blonr = fltarr(nz,nf,nl)
  openr,1,dir+'blonu.dat',/f77_unformatted
  print,' reading blon'
  readu,1,blonr

  blon = fltarr(nz,nf,nl+1)
  blon(*,*,0:nl-1) = blonr
  blon(*,*,nl)   = blon(*,*,0)

; magnetic latitude

  close,1
  blatr = fltarr(nz,nf,nl)
  openr,1,dir+'blatu.dat',/f77_unformatted
  print,' reading blat'
  readu,1,blatr

  blat = fltarr(nz,nf,nl+1)
  blat(*,*,0:nl-1) = blatr
  blat(*,*,nl)   = blat(*,*,0)

; altitude

  close,1
  baltr = fltarr(nz,nf,nl)
  openr,1,dir+'baltu.dat',/f77_unformatted
  print,' reading balt'
  readu,1,baltr

  balt = fltarr(nz,nf,nl+1)
;  balt(*,*,0:nl-1) = baltr
  balt(*,*,0:nl-1) = baltr - 6370.
  balt(*,*,nl)   = balt(*,*,0)

; coordinate space (?)

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

