; READ-U.PRO

; set the color table 
  
  device,true=24
  device,de=0
  loadct,15

; reads in selected data arrays from sami2-0.95a for unformatted data

; sami3 

  nz    = 101  
  nion  = 7    
  nf    = 130
  nl    = 96
  nt    = 50

  dir='/home/node19/col/spreadF/sami3_esf-1.00/tmp/'
  dir='/home/node19/col/spreadF/sami3_esf-1.00_pphi/'
  dir='/home/node19/col/spreadF/sami3_esf-1.00/run3/sami3_esf-1.00_zg/'
  dir='/home/ppdbw/col/sami3_esf-1.00/run1/'
  dir='/home/node19/col/spreadF/sami3_esf-1.00/run1/'
  dir='/home/node19/col/spreadF/sami3_esf-1.00_pphi/sami3_esf-1.00_pphi_zg/tanh_0vn/new_bc/'
  dir='/home/node19/col/spreadF/sami3_esf-1.00_pphi/sami3_esf-1.00_pphi_zg/tanh_0c/'
;  dir='/home/node19/col/spreadF/sami3_esf-1.00_pphi/sami3_esf-1.00_pphi_zg/'
;  dir='/home/node19/col/spreadF/sami3_esf-1.00_pphi/sami3_esf-1.00_pphi_zg_hall/'

  dir = '/raid2/col/spreadF_3D/sami3_esf-1.00_pphi/periodic/no_zonal_wind/no_hall/'
  dir='/home/node19/col/spreadF/sami3_esf-1.00_pphi/sami3_esf-1.00_pphi_p/tanh_0c/'
  dir='/home/node19/col/spreadF/sami3_esf-1.00_pphi/sami3_esf-1.00_pphi_p/tanh_0c/p_interpolate_spline/'

  dir='/home/node19/col/spreadF/sami3_esf-1.01p/'

 dir='/raid2/col/spreadF_3D/sami3_esf-1.03pz/run_020_020/'

; dir='//home/node24/col/spreadF_3D/sami3_esf-1.03pz/run_000h/'

 dir='/raid2/col/spreadF_3D/sami3_esf-1.00_pphi/periodic/zonal_wind/hall/'
 dir='/home/node27/col/esf_pphi/'

 dir='/home/node24/col/spreadF_3D/sami3_esf-1.10/'

  dt = .03333

  nx = 100  ; latitude
  ny = 100  ; altitude

  nnx = nl+1
  nny = nf-1

  nfp1 = nf + 1
; define directory of data files
; default is current directory

;  dir = './'


  t = 1
  tecc = 11
  irr0 = 11
  teccb = 0
  iden = 11
  idenn = 11
  iden1 = 11
  ivel = 11
  iti  = 0
  ite  = 11
  ivn  = 11
  idene= 1
  idene0 = 11
  idene0i = 0
  cont = 0
  magcoord = 1
  uu  = 1
  ip1 = 0
  coord = 11
  iphi  = 11
  rcm  = 11
  hic = 11
  hihc1 = 11
  sigmapp = 11
  sigmahp = 11
  phi1  = 1
  coef = 11

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

 if 1 eq 0 then begin

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
  denir    = fltarr(nz,nf,nl,7,nt)
  denitmp = fltarr(nz,nf,nl,7)
  openr,1,dir+'deniu.dat',/f77_unformatted
  print,' reading deniu '
  for i = 0,nt-1 do begin
    readu,1,denitmp
    denir(*,*,*,*,i) = denitmp
  endfor

  deni = fltarr(nz,nf,nl+1,7,nt)
  deni(*,*,0:nl-1,*,*) = denir
  deni(*,*,nl,*,*) = deni(*,*,0,*,*)

  dene = total(deni,4)

  endif

; ion density data

  if idenn eq 1 then begin

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
;  openr,1,dir+'deniu4.dat',/f77_unformatted
;  print,' reading deni4 '
;  for i = 0,nt-1 do begin
;    readu,1,denitmp
;    denir(*,*,*,i) = denitmp
;  endfor

;  deni4 = fltarr(nz,nf,nl+1,nt)
;  deni4(*,*,0:nl-1,*,*) = denir
;  deni4(*,*,nl,*,*) = deni4(*,*,0,*)

  close,1
  denir    = fltarr(nz,nf,nl,nt)
  denitmp = fltarr(nz,nf,nl)
  openr,1,dir+'deniu3.dat',/f77_unformatted
  print,' reading deni3 '
  for i = 0,nt-1 do begin
    readu,1,denitmp
    denir(*,*,*,i) = denitmp
  endfor

  deni3 = fltarr(nz,nf,nl+1,nt)
  deni3(*,*,0:nl-1,*,*) = denir
  deni3(*,*,nl,*,*) = deni3(*,*,0,*)


;  close,1
;  denir    = fltarr(nz,nf,nl,nt)
;  denitmp = fltarr(nz,nf,nl)
;  openr,1,dir+'deniu1.dat',/f77_unformatted
;  print,' reading deni1 '
;  for i = 0,nt-1 do begin
;    readu,1,denitmp
;    denir(*,*,*,i) = denitmp
;  endfor

;  deni1 = fltarr(nz,nf,nl+1,nt)
;  deni1(*,*,0:nl-1,*,*) = denir
;  deni1(*,*,nl,*,*) = deni1(*,*,0,*)

  close,1
  denir    = fltarr(nz,nf,nl,nt)
  denitmp = fltarr(nz,nf,nl)
  openr,1,dir+'deniu2.dat',/f77_unformatted
  print,' reading deni2 '
  for i = 0,nt-1 do begin
    readu,1,denitmp
    denir(*,*,*,i) = denitmp
  endfor

  deni2 = fltarr(nz,nf,nl+1,nt)
  deni2(*,*,0:nl-1,*,*) = denir
  deni2(*,*,nl,*,*) = deni2(*,*,0,*)


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

  endif

; ion temperature data

  if iti eq 1 then begin

  close,1
  tir    = fltarr(nz,nf,nl,nion,nt)
  titmp = fltarr(nz,nf,nl,nion)
  openr,1,dir+'tiu.dat',/f77_unformatted
  print,' reading ti '
  for i = 0,nt-1 do begin
    readu,1,titmp
    tir(*,*,*,*,i) = titmp
  endfor

  ti = fltarr(nz,nf,nl+1,nion,nt)
  ti(*,*,0:nl-1,*,*) = tir
  ti(*,*,nl,*,*) = ti(*,*,0,*,*)

  endif

; electron temperature data

  if ite eq 1 then begin

  close,1
  ter    = fltarr(nz,nf,nl,nt)
  tetmp = fltarr(nz,nf,nl)
  openr,1,dir+'teu.dat',/f77_unformatted
  print,' reading te '
  for i = 0,nt-1 do begin
    readu,1,tetmp
    ter(*,*,*,i) = tetmp
  endfor

  te = fltarr(nz,nf,nl+1,nt)
  te(*,*,0:nl-1,*) = ter
  te(*,*,nl,*) = te(*,*,0,*)

  endif

; neutral wind data

  if ivn eq 1 then begin

  close,1
  vnqr    = fltarr(nz,nf,nl,nt)
  vnqtmp = fltarr(nz,nf,nl)
  openr,1,dir+'vnqu.dat',/f77_unformatted
  print,' reading vnq '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,vnqtmp
    vnqr(*,*,*,i) = vnqtmp
  endfor

  vnq = fltarr(nz,nf,nl+1,nt)
  vnq(*,*,0:nl-1,*) = vnqr
  vnq(*,*,nl,*) = vnq(*,*,0,*)

  close,1
  vnpr    = fltarr(nz,nf,nl,nt)
  vnptmp = fltarr(nz,nf,nl)
  openr,1,dir+'vnpu.dat',/f77_unformatted
  print,' reading vnp '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,vnptmp
    vnpr(*,*,*,i) = vnptmp
  endfor



  vnp = fltarr(nz,nf,nl+1,nt)
  vnp(*,*,0:nl-1,*) = vnpr
  vnp(*,*,nl,*) = vnp(*,*,0,*)

  close,1
  vnphir    = fltarr(nz,nf,nl,nt)
  vnphitmp = fltarr(nz,nf,nl)
  openr,1,dir+'vnphiu.dat',/f77_unformatted
  print,' reading vnphi '
  for i = 0,nt-1 do begin
    print,i 
    readu,1,vnphitmp
    vnphir(*,*,*,i) = vnphitmp
  endfor



  vnphi = fltarr(nz,nf,nl+1,nt)
  vnphi(*,*,0:nl-1,*) = vnphir
  vnphi(*,*,nl,*) = vnphi(*,*,0,*)

  endif

; electron density data

  if idene eq 1 then begin
  a=long(1.)
  close,1
  dene    = fltarr(nz,nf,nl,nt)
  denetmp = fltarr(nz,nf,nl)
  openr,1,dir+'deneu.dat',/f77_unformatted
  print,' reading dene '
  for i = 0,nt-1 do begin
;    print,i 
    readu,1,denetmp
    dene(*,*,*,i) = denetmp
  endfor

;  delvar,dener

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

  endif

  if coord eq 1 then begin

; coordinate space

  close,1
  xr = fltarr(nz,nf,nl)
  openr,1,dir+'xu.dat',/f77_unformatted
  print,' reading x'
  readu,1,xr

  x = fltarr(nz,nf,nl+1)
  x(*,*,0:nl-1) = xr
  x(*,*,nl)   = x(*,*,0)

  close,1
  yr = fltarr(nz,nf,nl)
  openr,1,dir+'yu.dat',/f77_unformatted
  print,' reading y'
  readu,1,yr

  y = fltarr(nz,nf,nl+1)
  y(*,*,0:nl-1) = yr
  y(*,*,nl)   = y(*,*,0)

  close,1
  zr = fltarr(nz,nf,nl)
  openr,1,dir+'zu.dat',/f77_unformatted
  print,' reading z'
  readu,1,zr

  z = fltarr(nz,nf,nl+1)
  z(*,*,0:nl-1) = zr
  z(*,*,nl)   = z(*,*,0)

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
  u1=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u1u.dat',/f77_unformatted
  print,' reading u1'
  for i = 0,nt-1 do begin
    readu,1,etmp
    u1(*,*,*,i) = etmp
  endfor


;  delvar,u1r


; u2

;  close,1
;  u2r=fltarr(nz,nf,nl,nt)
;  etmp=fltarr(nz,nf,nl)
;  openr,1,dir+'u2u.dat',/f77_unformatted
;  print,' reading u2'
;  for i = 0,nt-1 do begin
;    readu,1,etmp
;    u2r(*,*,*,i) = etmp
;  endfor

;  u2   = fltarr(nz,nf,nl+1,nt)
;  u2(*,*,0:nl-1,*) = u2r
;  u2(*,*,nl,*)     = u2(*,*,0,*)

;  delvar,u2r

;  stop

; u3

  close,1
  u3=fltarr(nz,nf,nl,nt)
  etmp=fltarr(nz,nf,nl)
  openr,1,dir+'u3u.dat',/f77_unformatted
  print,' reading u3'
  for i = 0,nt-1 do begin
    readu,1,etmp
    u3(*,*,*,i) = etmp
  endfor

;  u3   = fltarr(nz,nf,nl+1,nt)
;  u3(*,*,0:nl-1,*) = u3r
;  u3(*,*,nl,*)     = u3(*,*,0,*)

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

; hipcp

  close,1
  hipcp = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hipcpu.dat',/f77_unformatted
  print,' reading hipcp'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hipcp(*,*,i) = hipcrt
  endfor

; hidphig

  close,1
  hidphig = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hidphigu.dat',/f77_unformatted
  print,' reading hidphig'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hidphig(*,*,i) = hipcrt
  endfor

; hipcphi

  close,1
  hipcphi = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hipcphiu.dat',/f77_unformatted
  print,' reading hipcphi'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hipcphi(*,*,i) = hipcrt
  endfor


  endif

  if hihc1 eq 1 then begin
 
; hihc

  close,1
  hihcm = fltarr(nf,nl,nt)
  hihcrt= fltarr(nf,nl)
  openr,1,dir+'hihcmu.dat',/f77_unformatted
  print,' reading hihcm'
  for i = 0,nt-1 do begin
    readu,1,hihcrt
    hihcm(*,*,i) = hihcrt
  endfor

;  hihcm = fltarr(nf,nl+1,nt)
;  hihcm(*,0:nl-1,*) = hihcr(*,*,*)
;  hihcm(*,nl,*)     = hihcr(*,0,*)

endif

  if 1 eq 11 then begin

; hidpv

  close,1

  hidpv = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hidpvu.dat',/f77_unformatted
  print,' reading hidpv'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hidpv(*,*,i) = hipcrt
  endfor

;  hidpv = fltarr(nf,nl+1,nt)
;  hidpv(*,0:nl-1,*) = hipcr(*,*,*)
;  hidpv(*,nl,*)     = hipcr(*,0,*)

  endif

  if 1 eq 0 then begin

; hidpg

  close,1
  hipcr = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hidpgu.dat',/f77_unformatted
  print,' reading hidpg'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hipcr(*,*,i) = hipcrt
  endfor

  hidpg = fltarr(nf,nl+1,nt)
  hidpg(*,0:nl-1,*) = hipcr(*,*,*)
  hidpg(*,nl,*)     = hipcr(*,0,*)

; hidphiv

  close,1
  hipcr = fltarr(nf,nl,nt)
  hipcrt= fltarr(nf,nl)
  openr,1,dir+'hidphivu.dat',/f77_unformatted
  print,' reading hidphiv'
  for i = 0,nt-1 do begin
    readu,1,hipcrt
    hipcr(*,*,i) = hipcrt
  endfor

  hidphiv = fltarr(nf,nl+1,nt)
  hidphiv(*,0:nl-1,*) = hipcr(*,*,*)
  hidphiv(*,nl,*)     = hipcr(*,0,*)


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
;    print,i 
    readu,1,sigmaptmp
    sigmapr(*,*,*,i) = sigmaptmp
  endfor

  sigmap = fltarr(nz,nf,nl+1,nt)
  sigmap(*,*,0:nl-1,*) = sigmapr
  sigmap(*,*,nl,*) = sigmap(*,*,0,*)

  close,1
  sigmapr    = fltarr(nz,nf,nl,nt)
  sigmaptmp = fltarr(nz,nf,nl)
  openr,1,dir+'sigmapicu.dat',/f77_unformatted
  print,' reading sigmapic '
  for i = 0,nt-1 do begin
;    print,i 
    readu,1,sigmaptmp
    sigmapr(*,*,*,i) = sigmaptmp
  endfor

  sigmapic = fltarr(nz,nf,nl+1,nt)
  sigmapic(*,*,0:nl-1,*) = sigmapr
  sigmapic(*,*,nl,*) = sigmapic(*,*,0,*)

;  delvar, sigmapr

  endif

  if sigmahp eq 1 then begin

  close,1
  sigmahr    = fltarr(nz,nf,nl,nt)
  sigmahtmp = fltarr(nz,nf,nl)
  openr,1,dir+'sigmahu.dat',/f77_unformatted
  print,' reading sigmah '
  for i = 0,nt-1 do begin
;    print,i 
    readu,1,sigmahtmp
    sigmahr(*,*,*,i) = sigmahtmp
  endfor

  sigmah = fltarr(nz,nf,nl+1,nt)
  sigmah(*,*,0:nl-1,*) = sigmahr
  sigmah(*,*,nl,*) = sigmah(*,*,0,*)

  close,1
  sigmahr    = fltarr(nz,nf,nl,nt)
  sigmahtmp = fltarr(nz,nf,nl)
  openr,1,dir+'sigmahicu.dat',/f77_unformatted
  print,' reading sigmahic '
  for i = 0,nt-1 do begin
;    print,i 
    readu,1,sigmahtmp
    sigmahr(*,*,*,i) = sigmahtmp
  endfor

  sigmahic = fltarr(nz,nf,nl+1,nt)
  sigmahic(*,*,0:nl-1,*) = sigmahr
  sigmahic(*,*,nl,*) = sigmahic(*,*,0,*)

;  delvar, sigmahr

  endif

  if phi1 eq 1 then begin

  close,1
  phi    = fltarr(nnx,nny,nt)
  phitmp = fltarr(nnx,nny)
  openr,1,dir+'phiu.dat',/f77_unformatted
  print,' reading phi'
  for i = 0,nt-1 do begin
;    print,i 
    readu,1,phitmp
    phi(*,*,i) = phitmp
  endfor


  endif

  if coef eq 1 then begin

  close,1
  cxxe    = fltarr(nnx,nny,nt)
  cyye    = fltarr(nnx,nny,nt)
  cxe     = fltarr(nnx,nny,nt)
  cye     = fltarr(nnx,nny,nt)
  rhse    = fltarr(nnx,nny,nt)
  phitmp1 = fltarr(nnx,nny)
  phitmp2 = fltarr(nnx,nny)
  phitmp3 = fltarr(nnx,nny)
  phitmp4 = fltarr(nnx,nny)
  phitmp5 = fltarr(nnx,nny)
  openr,1,dir+'c.dat',/f77_unformatted
  print,' reading c'
  for i = 0,nt-1 do begin
;    print,i 
    readu,1,phitmp1,phitmp2,phitmp3,phitmp4,phitmp5
    cxxe(*,*,i) = phitmp1
    cyye(*,*,i) = phitmp2
    cxe(*,*,i)  = phitmp3
    cye(*,*,i)  = phitmp4
    rhse(*,*,i) = phitmp5
  endfor


  endif


  close,1

    delvar,dener,u1r,u3r

  end
