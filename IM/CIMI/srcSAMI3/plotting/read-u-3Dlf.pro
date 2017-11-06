; READ-U-3DLF.PRO

; set the color table 
  
  device,true=24
  device,de=0
  device,retain=2
  loadct,41

; reads in selected data arrays from sami2-0.95a for unformatted data

; default:

;   time.dat
;   deniu.dat
;   glatu.dat
;   zaltu.dat

; uncomment other data arrays of interest

  nz    = 201
  nion  = 7
  nf    = 202
  nl    = 96
  nt   =  120

  nz    = 100
  nion  = 7
  nf    = 98
  nl    = 96
  nt   =  6


  nx = 200  ; latitude
  ny = 200  ; altitude

  nf2 = 126

  nnx = nl+1
  nny = nf-1

;  dir='/marsh/data2/col/ionosphere/esf/sami3_esf-1.50_vw/tilt_0/tilt_0_400/3D/'
  dir='/marsh/data2/col/ionosphere/esf/aveiro/code5_vn0_vexb0/'

  dt = .03333

  close,1

  t        = 1
  tecc     = 13
  idene    = 1
  phi1     = 14
  uu       = 14
  upsh     = 1
  ite      = 111
  iden     = 111
  hic      = 12
  sigmapp  = 11
  sigmahp  = 11
  ihidpv   = 12
  ivnp     = 11
  igeom    = 1
  ivel     = 144
  iti      = 111

  ifuv     = 11
  irr0     = 0
  teccb    = 0
  idenn    = 0
  iden1    = 0
  ivn      = 0
  idene0   = 321
  idene0B  = 13
  idene0i  = 0
  cont     = 0
  magcoord = 1
  ip1      = 11
  coord    = 1
  iphi     = 0
  rcm      = 0
  icurrent = 11
  icurrent0= 11
  ibfield  = 0



;  close,1
;  openr,1,dir+'bx0.dat',/f77_unformatted
;  bx0 = fltarr(nx,ny,nl)
;  readu,1,bx0

;  close,1
;  openr,1,dir+'by0.dat',/f77_unformatted
;  by0 = fltarr(nx,ny,nl)
;  readu,1,by0

;  close,1
;  openr,1,dir+'bz0.dat',/f77_unformatted
;  bz0 = fltarr(nx,ny,nl)
;  readu,1,bz0


  if ibfield eq 1 then begin
    close,1
    openr,1,dir+'bxu.dat',/f77_unformatted
    bx = fltarr(nz,nf,nl)
    readu,1,bx
    close,1
    openr,1,dir+'byu.dat',/f77_unformatted
    by = fltarr(nz,nf,nl)
    readu,1,by
    close,1
    openr,1,dir+'bzu.dat',/f77_unformatted
    bz = fltarr(nz,nf,nl)
    readu,1,bz
    close,1
    openr,1,dir+'bmst.dat',/f77_unformatted
    bms = fltarr(nz,nf,nl)
    readu,1,bms
  endif

  if icurrent eq 1 then begin

  close,1
  jpr    = fltarr(nz,nf,nl,nt)
  jtmp = fltarr(nz,nf,nl)
  openr,1,dir+'jpu.dat',/f77_unformatted
  print,' reading jp '
  for i = 0,nt-1 do begin
      print ,i
    readu,1,jtmp
    jpr(*,*,*,i) = jtmp
  endfor

  jp = fltarr(nz,nf,nl+1,nt)
  jp(*,*,0:nl-1,*) = jpr
  jp(*,*,nl,*) = jp(*,*,0,*)

  close,1
  jphir    = fltarr(nz,nf,nl,nt)
  jtmp = fltarr(nz,nf,nl)
  openr,1,dir+'jphiu.dat',/f77_unformatted
  print,' reading jphi '
  for i = 0,nt-1 do begin
    readu,1,jtmp
    jphir(*,*,*,i) = jtmp
  endfor

  jphi = fltarr(nz,nf,nl+1,nt)
  jphi(*,*,0:nl-1,*) = jphir
  jphi(*,*,nl,*) = jphi(*,*,0,*)

endif

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

;  phip = fltarr(0:nl-2,0:nf-2,nt)
  phip = reform(phi(0:nl-2,0:nf-2,*))

   close,1
  endif

; time data (formatted)

  tm = findgen(nt) * dt
  time = fltarr(4,nt)

  if irr0 eq 1 then begin

  close,1
  irrr = fltarr(nx,nl,nt)
  openr,1,dir+'irru.dat',/f77_unformatted
  print,' reading irr'
  readu,1,irrr

  irr = fltarr(nx,nl+1,nt)
  irr(*,0:nl-1,*) = irrr(*,*,*)
  irr(*,nl,*)     = irrr(*,0,*)

  endif

  close,1

  if ihidpv eq 1 then begin

  close,1
  hidpv = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hidpvu.dat',/f77_unformatted
  print,' reading hidpv'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hidpv(*,*,i) = hipcrt
  endfor

  close,1
  hidphiv = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hidphivu.dat',/f77_unformatted
  print,' reading hidphiv'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hidphiv(*,*,i) = hipcrt
  endfor

  endif

  if t eq 1 then begin

  close,1
  time = fltarr(4,nt)
  openr,1,dir+'time.dat'
  print,' reading time'
  readf,1,time

  endif

; tec,nmf2,hmf2,zalt0,glat0,glon0

  if tecc eq 1 then begin

  close,1
  glat0r = fltarr(nx,ny,nl)
  openr,1,dir+'glat0.dat',/f77_unformatted
  print,' reading glat0'
  readu,1,glat0r

  glat0 = fltarr(nx,ny,nl+1)
  glat0(*,*,0:nl-1 ) = glat0r(*,*,*)
  glat0(*,*,nl)      = glat0r(*,*,0)

  close,1
  glon0r = fltarr(nx,ny,nl)
  openr,1,dir+'glon0.dat',/f77_unformatted
  print,' reading glon0'
  readu,1,glon0r

  glon0 = fltarr(nx,ny,nl+1)
  glon0(*,*,0:nl-1 ) = glon0r(*,*,*)
  glon0(*,*,nl)      = glon0r(*,*,0)

  close,1
  zalt0r = fltarr(nx,ny,nl)
  openr,1,dir+'zalt0.dat',/f77_unformatted
  print,' reading zalt0'
  readu,1,zalt0r

  zalt0 = fltarr(nx,ny,nl+1)
  zalt0(*,*,0:nl-1 ) = zalt0r(*,*,*)
  zalt0(*,*,nl)      = zalt0r(*,*,0)

  close,1
  tecr = fltarr(nx,nl,nt)
  openr,1,dir+'tecu.dat',/f77_unformatted
  print,' reading tec'
  readu,1,tecr

  tec = fltarr(nx,nl+1,nt)
  tec(*,0:nl-1,*) = tecr(*,*,*)
  tec(*,nl,*)     = tecr(*,0,*)

  close,1
  nmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'nmf2u.dat',/f77_unformatted
  print,' reading nmf2'
  readu,1,nmf2r

  nmf2 = fltarr(nx,nl+1,nt)
  nmf2(*,0:nl-1,*) = nmf2r(*,*,*)
  nmf2(*,nl,*)     = nmf2r(*,0,*)

  close,1
  hmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'hmf2u.dat',/f77_unformatted
  print,' reading hmf2'
  readu,1,hmf2r

  hmf2 = fltarr(nx,nl+1,nt)
  hmf2(*,0:nl-1,*) = hmf2r(*,*,*)
  hmf2(*,nl,*)     = hmf2r(*,0,*)

endif

if ifuv eq 1 then begin

  close,1
  fuvr = fltarr(nx,nl,nt)
  openr,1,dir+'fuvu.dat',/f77_unformatted
  print,' reading fuv'
  readu,1,fuvr

  fuv = fltarr(nx,nl+1,nt)
  fuv(*,0:nl-1,*) = fuvr(*,*,*)
  fuv(*,nl,*)     = fuvr(*,0,*)

endif

; if 1 eq 0 then begin

  if idene0B eq 1 then begin

  close,1
  dener0    = fltarr(nx,ny,nl,nt)
  openr,1,dir+'dene0B.dat',/f77_unformatted
  print,' reading dene0B '
  readu,1,dener0

  dene0B = fltarr(nx,ny,nl+1,nt)
  dene0B(*,*,0:nl-1,*) = dener0
  dene0B(*,*,nl,*) = dene0B(*,*,0,*)

  close,1
  blat0r = fltarr(nx,ny,nl)
  openr,1,dir+'blat0.dat',/f77_unformatted
  print,' reading blat0'
  readu,1,blat0r

  blat0 = fltarr(nx,ny,nl+1)
  blat0(*,*,0:nl-1 ) = blat0r(*,*,*)
  blat0(*,*,nl)      = blat0r(*,*,0)

  close,1
  blon0r = fltarr(nx,ny,nl)
  openr,1,dir+'blon0.dat',/f77_unformatted
  print,' reading blon0'
  readu,1,blon0r

  blon0 = fltarr(nx,ny,nl+1)
  blon0(*,*,0:nl-1 ) = blon0r(*,*,*)
  blon0(*,*,nl)      = blon0r(*,*,0)

  close,1
  balt0r = fltarr(nx,ny,nl)
  openr,1,dir+'balt0.dat',/f77_unformatted
  print,' reading balt0'
  readu,1,balt0r

  balt0 = fltarr(nx,ny,nl+1)
  balt0(*,*,0:nl-1 ) = balt0r(*,*,*)
  balt0(*,*,nl)      = balt0r(*,*,0)

endif

  if teccb eq 1 then begin

  close,1
  blat0r = fltarr(nx,ny,nl)
  openr,1,dir+'blat0.dat',/f77_unformatted
  print,' reading blat0'
  readu,1,blat0r

  blat0 = fltarr(nx,ny,nl+1)
  blat0(*,*,0:nl-1 ) = blat0r(*,*,*)
  blat0(*,*,nl)      = blat0r(*,*,0)

  close,1
  blon0r = fltarr(nx,ny,nl)
  openr,1,dir+'blon0.dat',/f77_unformatted
  print,' reading blon0'
  readu,1,blon0r

  blon0 = fltarr(nx,ny,nl+1)
  blon0(*,*,0:nl-1 ) = blon0r(*,*,*)
  blon0(*,*,nl)      = blon0r(*,*,0)

  close,1
  balt0r = fltarr(nx,ny,nl)
  openr,1,dir+'balt0.dat',/f77_unformatted
  print,' reading balt0'
  readu,1,balt0r

  balt0 = fltarr(nx,ny,nl+1)
  balt0(*,*,0:nl-1 ) = balt0r(*,*,*)
  balt0(*,*,nl)      = balt0r(*,*,0)

  close,1
  tecr = fltarr(nx,nl,nt)
  openr,1,dir+'tecub.dat',/f77_unformatted
  print,' reading tecb'
  readu,1,tecr

  tecb = fltarr(nx,nl+1,nt)
  tecb(*,0:nl-1,*) = tecr(*,*,*)
  tecb(*,nl,*)     = tecr(*,0,*)

  close,1
  nmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'nmf2ub.dat',/f77_unformatted
  print,' reading nmf2b'
  readu,1,nmf2r

  nmf2b = fltarr(nx,nl+1,nt)
  nmf2b(*,0:nl-1,*) = nmf2r(*,*,*)
  nmf2b(*,nl,*)     = nmf2r(*,0,*)

  close,1
  hmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'hmf2ub.dat',/f77_unformatted
  print,' reading hmf2b'
  readu,1,hmf2r

  hmf2b = fltarr(nx,nl+1,nt)
  hmf2b(*,0:nl-1,*) = hmf2r(*,*,*)
  hmf2b(*,nl,*)     = hmf2r(*,0,*)

  endif

; ion density data

  if iden1 eq 1 then begin

  close,1
  denir    = fltarr(nz,nf,nl,nion,nt)
  denitmp = fltarr(nz,nf,nl,nion)
  openr,1,dir+'deniu.dat',/f77_unformatted
  print,' reading deniu '
  for i = 0,nt-1 do begin
    readu,1,denitmp
    denir(*,*,*,*,i) = denitmp
  endfor

  deni = fltarr(nz,nf,nl+1,nion,nt)
  deni(*,*,0:nl-1,*,*) = denir
  deni(*,*,nl,*,*) = deni(*,*,0,*,*)

  dene = total(deni,4)

  endif

; neutral density data

  if idenn eq 1 then begin

  close,1
  dennr    = fltarr(nz,nf,nl,nt)
  denntmp = fltarr(nz,nf,nl)
  openr,1,dir+'dennu2.dat',/f77_unformatted
  print,' reading denn1 '
  for i = 0,nt-1 do begin
    readu,1,denntmp
    dennr(*,*,*,i) = denntmp
  endfor

  denn2 = fltarr(nz,nf,nl+1,nt)
  denn2(*,*,0:nl-1,*,*) = dennr
  denn2(*,*,nl,*,*) = denn2(*,*,0,*)

  close,1
  dennr    = fltarr(nz,nf,nl,nt)
  denntmp = fltarr(nz,nf,nl)
  openr,1,dir+'dennu3.dat',/f77_unformatted
  print,' reading denn3 '
  for i = 0,nt-1 do begin
    readu,1,denntmp
    dennr(*,*,*,i) = denntmp
  endfor

  denn3 = fltarr(nz,nf,nl+1,nt)
  denn3(*,*,0:nl-1,*,*) = dennr
  denn3(*,*,nl,*,*) = denn3(*,*,0,*)

  close,1
  dennr    = fltarr(nz,nf,nl,nt)
  denntmp = fltarr(nz,nf,nl)
  openr,1,dir+'dennu4.dat',/f77_unformatted
  print,' reading denn4 '
  for i = 0,nt-1 do begin
    readu,1,denntmp
    dennr(*,*,*,i) = denntmp
  endfor

  denn4 = fltarr(nz,nf,nl+1,nt)
  denn4(*,*,0:nl-1,*,*) = dennr
  denn4(*,*,nl,*,*) = denn4(*,*,0,*)

  endif

; ion density data

  if iden eq 1 then begin

;  close,1
;  denir    = fltarr(nz,nf,nl,nt)
;  denitmp = fltarr(nz,nf,nl)
;  openr,1,dir+'deniu8.dat',/f77_unformatted
;  print,' reading deni8'
;  for i = 0,nt-1 do begin
;      print,i
;    readu,1,denitmp
;    denir(*,*,*,i) = denitmp
;  endfor

;  deni8 = fltarr(nz,nf,nl+1,nt)
;  deni8(*,*,0:nl-1,*,*) = denir
;  deni8(*,*,nl,*,*) = deni8(*,*,0,*)


;  stop

  close,1
  denir    = fltarr(nz,nf,nl,nt)
  denitmp = fltarr(nz,nf,nl) 
  openr,1,dir+'deniu1.dat',/f77_unformatted
  print,' reading deni1 '
  for i = 0,nt-1 do begin
     print,'deni1 - ',i
    readu,1,denitmp
    denir(*,*,*,i) = denitmp
  endfor

  deni1 = fltarr(nz,nf,nl+1,nt)
  deni1(*,*,0:nl-1,*,*) = denir
  deni1(*,*,nl,*,*) = deni1(*,*,0,*)

  close,1
  denir    = fltarr(nz,nf,nl,nt)
  denitmp = fltarr(nz,nf,nl)
  openr,1,dir+'deniu2.dat',/f77_unformatted
  print,' reading deni2 '
  for i = 0,nt-1 do begin
     print,'deni2 - ',i
    readu,1,denitmp
    denir(*,*,*,i) = denitmp
  endfor

  deni2 = fltarr(nz,nf,nl+1,nt)
  deni2(*,*,0:nl-1,*,*) = denir
  deni2(*,*,nl,*,*) = deni2(*,*,0,*)

  close,1
  denir    = fltarr(nz,nf,nl,nt)
  denitmp = fltarr(nz,nf,nl) 
  openr,1,dir+'deniu5.dat',/f77_unformatted
  print,' reading deni5 '
  for i = 0,nt-1 do begin
     print,'deni5 - ',i
    readu,1,denitmp
    denir(*,*,*,i) = denitmp
  endfor

  deni5 = fltarr(nz,nf,nl+1,nt)
  deni5(*,*,0:nl-1,*,*) = denir
  deni5(*,*,nl,*,*) = deni5(*,*,0,*)

  close,1
  denir    = fltarr(nz,nf,nl,nt)
  denitmp = fltarr(nz,nf,nl) 
  openr,1,dir+'deniu3.dat',/f77_unformatted
  print,' reading deni3 '
  for i = 0,nt-1 do begin
     print,'deni3 - ',i
    readu,1,denitmp
    denir(*,*,*,i) = denitmp
  endfor

  deni3 = fltarr(nz,nf,nl+1,nt)
  deni3(*,*,0:nl-1,*,*) = denir
  deni3(*,*,nl,*,*) = deni3(*,*,0,*)

  close,1
  denir    = fltarr(nz,nf,nl,nt)
  denitmp = fltarr(nz,nf,nl) 
  openr,1,dir+'deniu4.dat',/f77_unformatted
  print,' reading deni4 '
  for i = 0,nt-1 do begin
     print,'deni4 - ',i
    readu,1,denitmp
    denir(*,*,*,i) = denitmp
  endfor

  deni4 = fltarr(nz,nf,nl+1,nt)
  deni4(*,*,0:nl-1,*,*) = denir
  deni4(*,*,nl,*,*) = deni4(*,*,0,*)

  endif

; ion velocity data

  if ivel eq 1 then begin

  close,1
  vsir    = fltarr(nz,nf,nl,nt)
  vsitmp = fltarr(nz,nf,nl)
  openr,1,dir+'vsiu1.dat',/f77_unformatted
  print,' reading vsi1 '
  for i = 0,nt-1 do begin
    readu,1,vsitmp
    vsir(*,*,*,i) = vsitmp
  endfor

  vsi1 = fltarr(nz,nf,nl+1,nt)
  vsi1(*,*,0:nl-1,*) = vsir
  vsi1(*,*,nl,*) = vsi1(*,*,0,*)

  close,1
  vsir    = fltarr(nz,nf,nl,nt)
  vsitmp = fltarr(nz,nf,nl)
  openr,1,dir+'vsiu2.dat',/f77_unformatted
  print,' reading vsi2 '
  for i = 0,nt-1 do begin
    readu,1,vsitmp
    vsir(*,*,*,i) = vsitmp
  endfor

  vsi2 = fltarr(nz,nf,nl+1,nt)
  vsi2(*,*,0:nl-1,*) = vsir
  vsi2(*,*,nl,*) = vsi2(*,*,0,*)

  close,1
  vsir    = fltarr(nz,nf,nl,nt)
  vsitmp = fltarr(nz,nf,nl)
  openr,1,dir+'vsiu3.dat',/f77_unformatted
  print,' reading vsi3 '
  for i = 0,nt-1 do begin
    readu,1,vsitmp
    vsir(*,*,*,i) = vsitmp
  endfor

  vsi3 = fltarr(nz,nf,nl+1,nt)
  vsi3(*,*,0:nl-1,*) = vsir
  vsi3(*,*,nl,*) = vsi3(*,*,0,*)

  endif

; ion temperature data

  if iti eq 1 then begin

  close,1
  tir    = fltarr(nz,nf,nl,nt)
  titmp = fltarr(nz,nf,nl) 
  openr,1,dir+'tiu1.dat',/f77_unformatted
  print,' reading ti1 '
  for i = 0,nt-1 do begin
    readu,1,titmp
    tir(*,*,*,i) = titmp
  endfor

  ti1 = fltarr(nz,nf,nl+1,nt)
  ti1(*,*,0:nl-1,*,*) = tir
  ti1(*,*,nl,*,*) = ti1(*,*,0,*)

  close,1
  tir    = fltarr(nz,nf,nl,nt)
  titmp = fltarr(nz,nf,nl)
  openr,1,dir+'tiu2.dat',/f77_unformatted
  print,' reading ti2 '
  for i = 0,nt-1 do begin
    readu,1,titmp
    tir(*,*,*,i) = titmp
  endfor

  ti2 = fltarr(nz,nf,nl+1,nt)
  ti2(*,*,0:nl-1,*,*) = tir
  ti2(*,*,nl,*,*) = ti2(*,*,0,*)

  endif

; electron temperature data

  if ite eq 1 then begin

  close,1
  ter    = fltarr(nz,nf,nl,nt)
  tetmp = fltarr(nz,nf,nl)
  openr,1,dir+'teu.dat',/f77_unformatted
  print,' reading te '
  for i = 0,nt-1 do begin
     print ,i
    readu,1,tetmp
    ter(*,*,*,i) = tetmp
  endfor

  te = fltarr(nz,nf,nl+1,nt)
  te(*,*,0:nl-1,*) = ter
  te(*,*,nl,*) = te(*,*,0,*)

  endif

; neutral wind data magnetic

  if ivnp eq 1 then begin

  close,1
  vnp0    = fltarr(nz,nf,nl,nt)
  tetmp = fltarr(nz,nf,nl)
  openr,1,dir+'vnpu.dat',/f77_unformatted
  print,' reading vnpu.dat '
  for i = 0,nt-1 do begin
    readu,1,tetmp
    vnp0(*,*,*,i) = tetmp
  endfor

  vnp = fltarr(nz,nf,nl+1,nt)
  vnp(*,*,0:nl-1,*) = vnp0
  vnp(*,*,nl,*) = vnp0(*,*,0,*)

  close,1
  vnphi0    = fltarr(nz,nf,nl,nt)
  openr,1,dir+'vnphiu.dat',/f77_unformatted
  print,' reading vnphiu.dat '
  for i = 0,nt-1 do begin
    readu,1,tetmp
    vnphi0(*,*,*,i) = tetmp
  endfor

  vnphi = fltarr(nz,nf,nl+1,nt)
  vnphi(*,*,0:nl-1,*) = vnphi0
  vnphi(*,*,nl,*) = vnphi0(*,*,0,*)

  close,1
  vnq0    = fltarr(nz,nf,nl,nt)
  tetmp = fltarr(nz,nf,nl)
  openr,1,dir+'vnqu.dat',/f77_unformatted
  print,' reading vnqu.dat '
  for i = 0,nt-1 do begin
    readu,1,tetmp
    vnq0(*,*,*,i) = tetmp
  endfor

  vnq = fltarr(nz,nf,nl+1,nt)
  vnq(*,*,0:nl-1,*) = vnq0
  vnq(*,*,nl,*) = vnq0(*,*,0,*)

  endif



; neutral wind data

  if ivn eq 1 then begin

  close,1
  vnmerr    = fltarr(nz,nf,nl,nt)
  vtmp    = fltarr(nz,nf,nl)
  openr,1,dir+'vumer.dat',/f77_unformatted
  print,' reading vumer '
  for i = 0,nt-1 do begin
    readu,1,vtmp
    vnmerr(*,*,*,i) = vtmp
  endfor

  vnmer = fltarr(nz,nf,nl+1,nt)
  vnmer(*,*,0:nl-1,*) = vnmerr
  vnmer(*,*,nl,*) = vnmer(*,*,0,*)

  close,1
  vnzonr    = fltarr(nz,nf,nl,nt)
  vtmp    = fltarr(nz,nf,nl)
  openr,1,dir+'uuzon.dat',/f77_unformatted
  print,' reading uuzon '
  for i = 0,nt-1 do begin
    readu,1,vtmp
    vnzonr(*,*,*,i) = vtmp
  endfor

  vnzon = fltarr(nz,nf,nl+1,nt)
  vnzon(*,*,0:nl-1,*) = vnzonr
  vnzon(*,*,nl,*) = vnzon(*,*,0,*)

  endif

; electron density data

  if idene eq 1 then begin
  a=long(1.)
  close,1
  dener    = fltarr(nz,nf,nl,nt)
  denetmp = fltarr(nz,nf,nl)
  openr,1,dir+'deneu.dat',/f77_unformatted
  print,' reading dene '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,denetmp
    dener(*,*,*,i) = denetmp
  endfor



  dene = fltarr(nz,nf,nl+1,nt)
  dene(*,*,0:nl-1,*) = dener
  dene(*,*,nl,*) = dene(*,*,0,*)

;  delvar,dener

  endif

  if icurrent0 eq 1 then begin

  close,1
  j_merr0    = fltarr(nx,ny,nl,nt)
  openr,1,dir+'j_mer0.dat',/f77_unformatted
  print,' reading j_mer0 '
  readu,1,j_merr0

  j_mer0 = fltarr(nx,ny,nl+1,nt)
  j_mer0(*,*,0:nl-1,*) = j_merr0
  j_mer0(*,*,nl,*) = j_mer0(*,*,0,*)

  close,1
  j_zonr0    = fltarr(nx,ny,nl,nt)
  openr,1,dir+'j_zon0.dat',/f77_unformatted
  print,' reading j_zon0 '
  readu,1,j_zonr0

  j_zon0 = fltarr(nx,ny,nl+1,nt)
  j_zon0(*,*,0:nl-1,*) = j_zonr0
  j_zon0(*,*,nl,*) = j_zon0(*,*,0,*)

  close,1
  j_merr0    = fltarr(nx,nl,nt)
  openr,1,dir+'j_mer0_sq.dat',/f77_unformatted
  print,' reading j_mer0_sq '
  readu,1,j_merr0

  j_mer0_sq = fltarr(nx,nl+1,nt)
  j_mer0_sq(*,0:nl-1,*) = j_merr0
  j_mer0_sq(*,nl,*) = j_mer0_sq(*,0,*)

  close,1
  j_zonr0    = fltarr(nx,nl,nt)
  openr,1,dir+'j_zon0_sq.dat',/f77_unformatted
  print,' reading j_zon0_sq '
  readu,1,j_zonr0

  j_zon0_sq = fltarr(nx,nl+1,nt)
  j_zon0_sq(*,0:nl-1,*) = j_zonr0
  j_zon0_sq(*,nl,*) = j_zon0_sq(*,0,*)


  endif

  if idene0 eq 1 then begin

  close,1
  dener0    = fltarr(nx,ny,nl,nt)
  openr,1,dir+'dene0.dat',/f77_unformatted
  print,' reading dene0 '
  readu,1,dener0

  dene0 = fltarr(nx,ny,nl+1,nt)
  dene0(*,*,0:nl-1,*) = dener0
  dene0(*,*,nl,*) = dene0(*,*,0,*)

;  close,1
;  dener0    = fltarr(nx,ny,nl,nt)
;  openr,1,dir+'deni02.dat',/f77_unformatted
;  print,' reading deni02 '
;  readu,1,dener0

;  deni02 = fltarr(nx,ny,nl+1,nt)
;  deni02(*,*,0:nl-1,*) = dener0
;  deni02(*,*,nl,*) = deni02(*,*,0,*)

  close,1
  glat0r = fltarr(nx,ny,nl)
  openr,1,dir+'glat0.dat',/f77_unformatted
  print,' reading glat0'
  readu,1,glat0r

  glat0 = fltarr(nx,ny,nl+1)
  glat0(*,*,0:nl-1 ) = glat0r(*,*,*)
  glat0(*,*,nl)      = glat0r(*,*,0)

  close,1
  zalt0r = fltarr(nx,ny,nl)
  openr,1,dir+'zalt0.dat',/f77_unformatted
  print,' reading zalt0'
  readu,1,zalt0r

  zalt0 = fltarr(nx,ny,nl+1)
  zalt0(*,*,0:nl-1 ) = zalt0r(*,*,*)
  zalt0(*,*,nl)      = zalt0r(*,*,0)


; geophyscial longitude

  close,1
  glonr0 = fltarr(nx,ny,nl)
  openr,1,dir+'glon0.dat',/f77_unformatted
  print,' reading glon0'
  readu,1,glonr0

  glon0 = fltarr(nx,ny,nl+1)
  glon0(*,*,0:nl-1) = glonr0
  glon0(*,*,nl)   = glon0(*,*,0)

; geophysical latitude

;  close,1
;  blatr0 = fltarr(nx,ny,nl)
;  openr,1,dir+'blat0.dat',/f77_unformatted
;  print,' reading blat0'
;  readu,1,blatr0

;  blat0 = fltarr(nx,ny,nl+1)
;  blat0(*,*,0:nl-1) = blatr0
;  blat0(*,*,nl)   = blat0(*,*,0)

; altitude

  close,1
  zaltr0 = fltarr(nx,ny,nl)
  openr,1,dir+'zalt0.dat',/f77_unformatted
  print,' reading zalt0'
  readu,1,zaltr0

  zalt0 = fltarr(nx,ny,nl+1)
  zalt0(*,*,0:nl-1) = zaltr0
  zalt0(*,*,nl)   = zalt0(*,*,0)

  endif

  if idene0i eq 1 then begin

  close,1
  dener0i    = fltarr(nx,ny,nl,nt)
  openr,1,dir+'dene0i.dat',/f77_unformatted
  print,' reading dene0i '
  readu,1,dener0i

  dene0i = fltarr(nx,ny,nl+1,nt)
  dene0i(*,*,0:nl-1,*) = dener0i
  dene0i(*,*,nl,*) = dene0i(*,*,0,*)

; geophyscial longitude

  close,1
  glonr0i = fltarr(nx,ny,nl)
  openr,1,dir+'glon0i.dat',/f77_unformatted
  print,' reading glon0i'
  readu,1,glonr0i

  glon0i = fltarr(nx,ny,nl+1)
  glon0i(*,*,0:nl-1) = glonr0i
  glon0i(*,*,nl)   = glon0i(*,*,0)

; geophysical latitude

  close,1
  glatr0i = fltarr(nx,ny,nl)
  openr,1,dir+'glat0i.dat',/f77_unformatted
  print,' reading glat0i'
  readu,1,glatr0i

  glat0i = fltarr(nx,ny,nl+1)
  glat0i(*,*,0:nl-1) = glatr0i
  glat0i(*,*,nl)   = glat0i(*,*,0)

; altitude

  close,1
  zaltr0i = fltarr(nx,ny,nl)
  openr,1,dir+'zalt0i.dat',/f77_unformatted
  print,' reading zalt0i'
  readu,1,zaltr0i

  zalt0i = fltarr(nx,ny,nl+1)
  zalt0i(*,*,0:nl-1) = zaltr0i
  zalt0i(*,*,nl)   = zalt0i(*,*,0)

  endif

; geophyscial longitude

  close,1
  glonr = fltarr(nz,nf,nl)
  openr,1,dir+'glonu.dat',/f77_unformatted
  print,' reading glon'
  readu,1,glonr

  glon = fltarr(nz,nf,nl+1)
  glon(*,*,0:nl-1) = glonr
  glon(*,*,nl)   = glon(*,*,0)

; geophysical latitude

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

if magcoord eq 1 then begin

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
  balt(*,*,0:nl-1) = baltr
  balt(*,*,nl)   = balt(*,*,0)

  lonp = reform(blon(nz/2,0,0:nl-2))
  altp = reform(balt(nz/2,0:nf-2,0))-6370.

  endif

  if coord  eq 1 then begin

; coordinate space

  close,1
  xr = fltarr(nz,nf,nl)
  openr,1,dir+'xsu.dat',/f77_unformatted
  print,' reading x'
  readu,1,xr

  x = fltarr(nz,nf,nl+1)
  x(*,*,0:nl-1) = xr
  x(*,*,nl)   = x(*,*,0)

  close,1
  yr = fltarr(nz,nf,nl)
  openr,1,dir+'ysu.dat',/f77_unformatted
  print,' reading y'
  readu,1,yr

  y = fltarr(nz,nf,nl+1)
  y(*,*,0:nl-1) = yr
  y(*,*,nl)   = y(*,*,0)

  close,1
  zr = fltarr(nz,nf,nl)
  openr,1,dir+'zsu.dat',/f77_unformatted
  print,' reading z'
  readu,1,zr

  z = fltarr(nz,nf,nl+1)
  z(*,*,0:nl-1) = zr
  z(*,*,nl)   = z(*,*,0)

endif

  if upsh eq 1 then begin

; u1p

  close,1
  u1pr=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u1pu.dat',/f77_unformatted
  print,' reading u1p'
  for i = 0,nt-1 do begin
    readu,1,etmp
    u1pr(*,*,*,i) = etmp
  endfor

  u1p   = fltarr(nz,nf,nl+1,nt)
  u1p(*,*,0:nl-1,*) = u1pr
  u1p(*,*,nl,*)     = u1p(*,*,0,*)

; u2s

;  close,1
;  u2sr=fltarr(nz,nf,nl,nt)
;  etmp=fltarr(nz,nf,nl)
;  openr,1,dir+'u2su.dat',/f77_unformatted
;  print,' reading u2s'
;  for i = 0,nt-1 do begin
;    readu,1,etmp
;    u2sr(*,*,*,i) = etmp
;  endfor

;  u2s   = fltarr(nz,nf,nl+1,nt)
;  u2s(*,*,0:nl-1,*) = u2sr
;  u2s(*,*,nl,*)     = u2s(*,*,0,*)

; u3h

  close,1
  u3hr=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u3hu.dat',/f77_unformatted
  print,' reading u3h'
  for i = 0,nt-1 do begin
    readu,1,etmp
    u3hr(*,*,*,i) = etmp
  endfor

  u3h   = fltarr(nz,nf,nl+1,nt)
  u3h(*,*,0:nl-1,*) = u3hr
  u3h(*,*,nl,*)     = u3h(*,*,0,*)

  endif


; u

  if uu eq 1 then begin

; u4

;  close,1
;  u4r=fltarr(nz,nf,nl,nt)
;  etmp=fltarr(nz,nf,nl)
;  openr,1,dir+'u4u.dat',/f77_unformatted
;  print,' reading u4'
;  for i = 0,nt-1 do begin
;    readu,1,etmp
;    u4r(*,*,*,i) = etmp
;  endfor

;  u4   = fltarr(nz,nf,nl+1,nt)
;  u4(*,*,0:nl-1,*) = u4r
;  u4(*,*,nl,*)     = u4(*,*,0,*)

;  stop

; u1

  close,1
  u1r=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u1u.dat',/f77_unformatted
  print,' reading u1'
  for i = 0,nt-1 do begin
    readu,1,etmp
    u1r(*,*,*,i) = etmp
  endfor

  u1   = fltarr(nz,nf,nl+1,nt)
  u1(*,*,0:nl-1,*) = u1r
  u1(*,*,nl,*)     = u1(*,*,0,*)


; if 4 eq 1 then begin

;  delvar,u1r


; u2

  close,1
  u2r=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u2u.dat',/f77_unformatted
  print,' reading u2'
  for i = 0,nt-1 do begin
    readu,1,etmp
    u2r(*,*,*,i) = etmp
  endfor

  u2   = fltarr(nz,nf,nl+1,nt)
  u2(*,*,0:nl-1,*) = u2r
  u2(*,*,nl,*)     = u2(*,*,0,*)

;  delvar,u2r

;  stop

; u3

  close,1
  u3r=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u3u.dat',/f77_unformatted
  print,' reading u3'
  for i = 0,nt-1 do begin
    readu,1,etmp
    u3r(*,*,*,i) = etmp
  endfor

  u3   = fltarr(nz,nf,nl+1,nt)
  u3(*,*,0:nl-1,*) = u3r
  u3(*,*,nl,*)     = u3(*,*,0,*)

; u4

;  close,1
;  u4r=fltarr(nz,nf,nl,nt)
;  etmp=fltarr(nz,nf,nl)
;  openr,1,dir+'u4u.dat',/f77_unformatted
;  print,' reading u4'
;  for i = 0,nt-1 do begin
;    readu,1,etmp
;    u4r(*,*,*,i) = etmp
;  endfor

;  u4   = fltarr(nz,nf,nl+1,nt)
;  u4(*,*,0:nl-1,*) = u4r
;  u4(*,*,nl,*)     = u4(*,*,0,*)

; u5

;  close,1
;  u5r=fltarr(nz,nf,nl,nt)
;  etmp=fltarr(nz,nf,nl)
;  openr,1,dir+'u5u.dat',/f77_unformatted
;  print,' reading u5'
;  for i = 0,nt-1 do begin
;    readu,1,etmp
;    u5r(*,*,*,i) = etmp
;  endfor

;  u5   = fltarr(nz,nf,nl+1,nt)
;  u5(*,*,0:nl-1,*) = u5r
;  u5(*,*,nl,*)     = u5(*,*,0,*)

  endif

; phi

  if iphi eq 1 then begin

  close,1
  phi1r=fltarr(nf,nl,nt)
  phi1rt=fltarr(nf,nl)
  openr,1,dir+'hipcu.dat',/f77_unformatted
  print,' reading hipcu'
  for i = 0,nt-1 do begin
    readu,1,phi1rt
    phi1r(*,*,i) = phi1rt
  endfor

  phi1 = fltarr(nf,nl+1,nt)
  phi1(*,0:nl-1,*) = phi1r
  phi1(*,nl,*)     = phi1(*,0,*)

  close,1
  phi2r=fltarr(nf,nl,nt)
  phi2rt=fltarr(nf,nl)
  openr,1,dir+'hihcu.dat',/f77_unformatted
  print,' reading hihcu'
  for i = 0,nt-1 do begin
    readu,1,phi2rt
    phi2r(*,*,i) = phi2rt
  endfor

  phi2 = fltarr(nf,nl+1,nt)
  phi2(*,0:nl-1,*) = phi2r
  phi2(*,nl,*)     = phi2(*,0,*)

  endif

; p1

  if ip1 eq 1 then begin

  close,1
  p1r=fltarr(nfp1,nl,nt)
  p1tmp=fltarr(nfp1,nl)
  openr,1,dir+'p1u.dat',/f77_unformatted
  print,' reading p1u'
  for i = 0,nt-1 do begin
    readu,1,p1tmp
    p1r(*,*,i) = p1tmp
  endfor

  p1   = fltarr(nfp1,nl+1,nt)
  p1(*,0:nl-1,*) = p1r
  p1(*,nl,*)     = p1(*,0,*)

  close,1
  p2r=fltarr(nf,nl,nt)
  p2tmp=fltarr(nf,nl)
  openr,1,dir+'p2u.dat',/f77_unformatted
  print,' reading p2u'
  for i = 0,nt-1 do begin
    readu,1,p2tmp
    p2r(*,*,i) = p2tmp
  endfor

  p2   = fltarr(nf,nl+1,nt)
  p2(*,0:nl-1,*) = p2r
  p2(*,nl,*)     = p2(*,0,*)

  close,1
  p3r=fltarr(nf,nl,nt)
  p3tmp=fltarr(nf,nl)
  openr,1,dir+'p3u.dat',/f77_unformatted
  print,' reading p3u'
  for i = 0,nt-1 do begin
    readu,1,p3tmp
    p3r(*,*,i) = p3tmp
  endfor

  p3   = fltarr(nf,nl+1,nt)
  p3(*,0:nl-1,*) = p3r
  p3(*,nl,*)     = p3(*,0,*)

  endif

; ep

  if cont eq 1 then begin

  close,1
  epr=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'ep1u.dat',/f77_unformatted
  print,' reading ep'
  for i = 0,nt-1 do begin
    readu,1,etmp
    epr(*,*,*,i) = etmp
  endfor

  ep   = fltarr(nz,nf,nl+1,nt)
  ep(*,*,0:nl-1,*) = epr
  ep(*,*,nl,*)     = ep(*,*,0,*)

; eh

  close,1
  ehr=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'eh1u.dat',/f77_unformatted
  print,' reading eh'
  for i = 0,nt-1 do begin
    readu,1,etmp
    ehr(*,*,*,i) = etmp
  endfor

  eh   = fltarr(nz,nf,nl+1,nt)
  eh(*,*,0:nl-1,*) = ehr
  eh(*,*,nl,*)     = eh(*,*,0,*)

  endif

  if hic eq 1 then begin

; hipc

  close,1
  hipcr = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hipcpu.dat',/f77_unformatted
  print,' reading hipcp'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hipcr(*,*,i) = hipcrt
  endfor

  hipcp = fltarr(nf,nl+1,nt)
  hipcp(*,0:nl-1,*) = hipcr(*,*,*)
  hipcp(*,nl,*)     = hipcr(*,0,*)

; hihc


  close,1
  hihcr = fltarr(nf,nl,nt)
  hihcrt= fltarr(nf,nl)
  openr,1,dir+'hihcmu.dat',/f77_unformatted
  print,' reading hihcm'
  for i = 0,nt-1 do begin
    readu,1,hihcrt
    hihcr(*,*,i) = hihcrt
  endfor

  hihcm = fltarr(nf,nl+1,nt)
  hihcm(*,0:nl-1,*) = hihcr(*,*,*)
  hihcm(*,nl,*)     = hihcr(*,0,*)

  endif

  if rcm eq 1 then begin


; hipc_rcm

  close,1
  hipc_rcm=fltarr(155,99,nt)
  hipc_rcm_tmp=fltarr(155,99)
  openr,1,dir+'hipcu_rcm.dat',/f77_unformatted
  print,' reading hipc_rcm'
  for i = 0,nt-1 do begin
    readu,1,hipc_rcm_tmp
    hipc_rcm(*,*,i) = hipc_rcm_tmp
  endfor

; hihc_rcm

  close,1
  hihc_rcm=fltarr(155,99,nt)
  hihc_rcm_tmp=fltarr(155,99)
  openr,1,dir+'hihcu_rcm.dat',/f77_unformatted
  print,' reading hihc_rcm'
  for i = 0,nt-1 do begin
    readu,1,hihc_rcm_tmp
    hihc_rcm(*,*,i) = hihc_rcm_tmp
  endfor

; glat_rcm

  close,1
  glat_rcm = fltarr(155)
  openr,1,dir+'glat_rcmu.dat',/f77_unformatted
  print,' reading glat_rcm'
  readu,1,glat_rcm

; glon_rcm

  close,1
  glon_rcm = fltarr(99)
  openr,1,dir+'glon_rcmu.dat',/f77_unformatted
  print,' reading glonx_rcm'
  readu,1,glon_rcm

; ion temperature

;  close,1
;  ti   = fltarr(nz,nf,nion,nt)
;  openr,1,dir+'tiu.dat',/f77_unformatted
;  print,' reading ti'
;  readu,1,ti

; electron temperature

;  close,1
;  te   = fltarr(nz,nf,nt)
;  openr,1,dir+'teu.dat',/f77_unformatted
;  print,' reading te'
;  readu,1,te

; ion velocity along the geomagnetic field

;  close,1
;  vsi  = fltarr(nz,nf,nion,nt)
;  openr,1,dir+'vsiu.dat',/f77_unformatted
;  print,' reading vsi'
;  readu,1,vsi

; neutral density

;  close,1
;  denn  = fltarr(nz,nf,nion,nt)
;  openr,1,dir+'denn.dat',/f77_unformatted
;  print,' reading denn'
;  readu,1,denn

; assorted diagnostic arrays

;  close,1
;  u1    = fltarr(nz,nf,nt)
;  openr,1,dir+'u1u.dat',/f77_unformatted
;  print,' reading u1'
;  readu,1,u1

;  close,1
;  u2    = fltarr(nz,nf,nl,nt)
;  openr,1,dir+'u2u.dat',/f77_unformatted
;  print,' reading u2'
;  readu,1,u2

;  close,1
;  u3  = fltarr(nz,nf,nt)
;  openr,1,dir+'u3u.dat',/f77_unformatted
;  print,' reading u3'
;  readu,1,u3

;  close,1
;  u4  = fltarr(nz,nf,nt)
;  openr,1,dir+'u4u.dat',/f77_unformatted
;  print,' reading u4'
;  readu,1,u4

;  close,1
;  u5  = fltarr(nz,nf,nt)
;  openr,1,dir+'u5u.dat',/f77_unformatted
;  print,' reading u5'
;  readu,1,u5

;  close,1
;  t1  = fltarr(nz,nf,nion,nt)
;  openr,1,dir+'t1u.dat',/f77_unformatted
;  print,' reading t1'
;  readu,1,t1

;  close,1
;  t2 = fltarr(nz,nf,nion,nt)
;  openr,1,dir+'t2u.dat',/f77_unformatted
;  print,' reading t2'
;  readu,1,t2

;  close,1
;  t3 = fltarr(nz,nf,nion,nt)
;  openr,1,dir+'t3u.dat',/f77_unformatted
;  print,' reading t3'
;  readu,1,t3

   endif

; sigma

  if sigmapp eq 1 then begin

  close,1
  sigmapr    = fltarr(nz,nf,nl,nt)
  sigmaptmp = fltarr(nz,nf,nl)
  openr,1,dir+'sigmapu.dat',/f77_unformatted
  print,' reading sigmap '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,sigmaptmp
    sigmapr(*,*,*,i) = sigmaptmp
  endfor

  sigmap = fltarr(nz,nf,nl+1,nt)
  sigmap(*,*,0:nl-1,*) = sigmapr
  sigmap(*,*,nl,*) = sigmap(*,*,0,*)

;  close,1
;  sigmapr    = fltarr(nz,nf,nl,nt)
;  sigmaptmp = fltarr(nz,nf,nl)
;  openr,1,dir+'sigmapBacu.dat',/f77_unformatted
;  print,' reading sigmapBa '
;  for i = 0,nt-1 do begin
;    print,i 
;    readu,1,sigmaptmp
;    sigmapr(*,*,*,i) = sigmaptmp
;  endfor

;  sigmapBa = fltarr(nz,nf,nl+1,nt)
;  sigmapBa(*,*,0:nl-1,*) = sigmapr
;  sigmapBa(*,*,nl,*) = sigmapBa(*,*,0,*)

;  delvar, sigmapr

  endif

  if sigmahp eq 1 then begin

  close,1
  sigmahr    = fltarr(nz,nf,nl,nt)
  sigmahtmp = fltarr(nz,nf,nl)
  openr,1,dir+'sigmahu.dat',/f77_unformatted
  print,' reading sigmah '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,sigmahtmp
    sigmahr(*,*,*,i) = sigmahtmp
  endfor

  sigmah = fltarr(nz,nf,nl+1,nt)
  sigmah(*,*,0:nl-1,*) = sigmahr
  sigmah(*,*,nl,*) = sigmah(*,*,0,*)

;  delvar, sigmahr

  endif

; igeom

  if igeom eq 1 then begin

    close,1
    vpnxr = fltarr(nz+1,nf+1,nl)
    openr,1,dir+'vpnxu.dat',/f77_unformatted
    print,' reading vpnxu.dat'
    readu,1,vpnxr
    vpnx1 = fltarr(nz,nf,nl)
    vpnx1 = vpnxr(0:nz-1,0:nf-1,*)
    vpnx = fltarr(nz,nf,nl+1)
    vpnx(*,*,0:nl-1) = vpnx1
    vpnx(*,*,nl)     = vpnx(*,*,0)

    close,1
    vpnyr = fltarr(nz+1,nf+1,nl)
    openr,1,dir+'vpnyu.dat',/f77_unformatted
    print,' reading vpnyu.dat'
    readu,1,vpnyr
    vpny1 = fltarr(nz,nf,nl)
    vpny1 = vpnyr(0:nz-1,0:nf-1,*)
    vpny = fltarr(nz,nf,nl+1)
    vpny(*,*,0:nl-1) = vpny1
    vpny(*,*,nl)     = vpny(*,*,0)

    close,1
    vpnzr = fltarr(nz+1,nf+1,nl)
    openr,1,dir+'vpnzu.dat',/f77_unformatted
    print,' reading vpnzu.dat'
    readu,1,vpnzr
    vpnz1 = fltarr(nz,nf,nl)
    vpnz1 = vpnzr(0:nz-1,0:nf-1,*)
    vpnz = fltarr(nz,nf,nl+1)
    vpnz(*,*,0:nl-1) = vpnz1
    vpnz(*,*,nl)     = vpnz(*,*,0)

    close,1
    vhnxr = fltarr(nz+1,nf+1,nl)
    openr,1,dir+'vhnxu.dat',/f77_unformatted
    print,' reading vhnxu.dat'
    readu,1,vhnxr
    vhnx1 = fltarr(nz,nf,nl)
    vhnx1 = vhnxr(0:nz-1,0:nf-1,*)
    vhnx = fltarr(nz,nf,nl+1)
    vhnx(*,*,0:nl-1) = vhnx1
    vhnx(*,*,nl)     = vhnx(*,*,0)

    close,1
    vhnyr = fltarr(nz+1,nf+1,nl)
    openr,1,dir+'vhnyu.dat',/f77_unformatted
    print,' reading vhnyu.dat'
    readu,1,vhnyr
    vhny1 = fltarr(nz,nf,nl)
    vhny1 = vhnyr(0:nz-1,0:nf-1,*)
    vhny = fltarr(nz,nf,nl+1)
    vhny(*,*,0:nl-1) = vhny1
    vhny(*,*,nl)     = vhny(*,*,0)

    close,1
    vhnzr = fltarr(nz+1,nf+1,nl)
    openr,1,dir+'vhnzu.dat',/f77_unformatted
    print,' reading vhnzu.dat'
    readu,1,vhnzr
    vhnz1 = fltarr(nz,nf,nl)
    vhnz1 = vhnzr(0:nz-1,0:nf-1,*)
    vhnz = fltarr(nz,nf,nl+1)
    vhnz(*,*,0:nl-1) = vhnz1
    vhnz(*,*,nl)     = vhnz(*,*,0)

    close,1
    bdirsxr = fltarr(nz+1,nf,nl)
    openr,1,dir+'bdirsxu.dat',/f77_unformatted
    print,' reading bdirsxu.dat'
    readu,1,bdirsxr
    bdirsx1 = fltarr(nz,nf,nl)
    bdirsx1 = bdirsxr(0:nz-1,*,*)
    bdirsx  = fltarr(nz,nf,nl+1)
    bdirsx(*,*,0:nl-1) = bdirsx1
    bdirsx(*,*,nl)     = bdirsx(*,*,0) 

    close,1
    bdirsyr = fltarr(nz+1,nf,nl)
    openr,1,dir+'bdirsyu.dat',/f77_unformatted
    print,' reading bdirsyu.dat'
    readu,1,bdirsyr
    bdirsy1 = fltarr(nz,nf,nl)
    bdirsy1 = bdirsyr(0:nz-1,*,*)
    bdirsy  = fltarr(nz,nf,nl+1)
    bdirsy(*,*,0:nl-1) = bdirsy1
    bdirsy(*,*,nl)     = bdirsy(*,*,0) 

    close,1
    bdirszr = fltarr(nz+1,nf,nl)
    openr,1,dir+'bdirszu.dat',/f77_unformatted
    print,' reading bdirszu.dat'
    readu,1,bdirszr
    bdirsz1 = fltarr(nz,nf,nl)
    bdirsz1 = bdirszr(0:nz-1,*,*)
    bdirsz  = fltarr(nz,nf,nl+1)
    bdirsz(*,*,0:nl-1) = bdirsz1
    bdirsz(*,*,nl)     = bdirsz(*,*,0) 

  endif

  end





