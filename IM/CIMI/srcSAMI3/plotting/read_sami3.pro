
; READ-U-3DLF.PRO

; set the color table 
  
  device,true=24
  device,de=0
  device,retain=2
  loadct,39

; set array sizes

  nz    = 120
  nf    = 160
  nl    =  90
  nion  =   7
  nt    = 191

  nnx   = nl + 1
  nny   = nf - 1

  dir = './'
;  dir = '/data2/col/ionosphere/sami3/sami3_mpi-1.70/euvac_hwm93/test_run_day91/'

  dt = .25

  close,1

  t        = 1
  idene    = 1
  phi1     = 1
  upsh     = 1

  if phi1 eq 1 then begin

  close,1
  phi    = fltarr(nnx,nny,nt)
  phitmp = fltarr(nnx,nny)
  openr,1,dir+'phiu.dat',/f77_unformatted
  print,' reading phi'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,phitmp
    phi(*,*,i) = phitmp
  endfor

  endif

; time data (formatted)

  tm   = findgen(nt) * dt
  time = fltarr(4,nt)

  if t eq 1 then begin

  close,1
  time = fltarr(4,nt)
  openr,1,dir+'time.dat'
  print,' reading time'
  readf,1,time

  endif

; electron density data

  if idene eq 1 then begin
  a=long(1.)
  close,1
  dene    = fltarr(nz,nf,nl,nt)
  etmp = fltarr(nz,nf,nl)
  openr,1,dir+'deneu.dat',/f77_unformatted
  print,' reading dene '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,etmp
    dene(*,*,*,i) = etmp
  endfor

; geophyscial longitude

  close,1
  glon = fltarr(nz,nf,nl)
  openr,1,dir+'glonu.dat',/f77_unformatted
  print,' reading glon'
  readu,1,glon

; geophysical latitude

  close,1
  glat = fltarr(nz,nf,nl)
  openr,1,dir+'glatu.dat',/f77_unformatted
  print,' reading glat'
  readu,1,glat

; altitude

  close,1
  zalt = fltarr(nz,nf,nl)
  openr,1,dir+'zaltu.dat',/f77_unformatted
  print,' reading zalt'
  readu,1,zalt

  endif

  if upsh eq 1 then begin

; u1p

  close,1
  u1p=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u1pu.dat',/f77_unformatted
  print,' reading u1p'
  for i = 0,nt-1 do begin
    print,i
    readu,1,etmp
    u1p(*,*,*,i) = etmp
  endfor

; u3h

  close,1
  u3h=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u3hu.dat',/f77_unformatted
  print,' reading u3h'
  for i = 0,nt-1 do begin
    print,i
    readu,1,etmp
    u3h(*,*,*,i) = etmp
  endfor

  endif

  end




