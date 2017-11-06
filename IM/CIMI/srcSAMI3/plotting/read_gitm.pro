
; set the color table 
  
  device,true=24
  device,de=0
  device,retain=2
  loadct,41

  j0 = 0L
  j1 = 0L
  j2 = 0L
  j3 = 0L

; array sizes

  jlon  = 72
  jlat  = 72
  jalt  = 49
  jtime = 97

  glatt  = fltarr(jlat)
  glont  = fltarr(jlon)
  zt     = fltarr(jlon,jlat,jalt,jtime)

  unt   = fltarr(jlon,jlat,jalt,jtime)
  vntt  = fltarr(jlon,jlat,jalt,jtime)

  ttn   = fltarr(jlon,jlat,jalt,jtime)
  tn2   = fltarr(jlon,jlat,jalt,jtime)
  to2   = fltarr(jlon,jlat,jalt,jtime)
  to1   = fltarr(jlon,jlat,jalt,jtime)
  tn4s  = fltarr(jlon,jlat,jalt,jtime)
  tno   = fltarr(jlon,jlat,jalt,jtime)
  thyd  = fltarr(jlon,jlat,jalt,jtime)

  dir   = '/home/jack/huba/ionosphere/sami3gitm/gitm_data/'

   print,jlon,jlat,jalt,jtime

  close,1
  openr,1,dir+'GITM_IDL.inp',/f77_unformatted
;  readu,1,j0,j1,j2,j3
  readu,1,j0,j1,j2,j3,glatt,glont,zt,$
          unt,vntt,ttn,tn2,to2,to1,tn4s,tno,thyd
   close,1

   print,j0,j1,j2,j3

  end

