; MAP8-PHI.PRO

; plots phi

  set_plot,'x'
 
; set the window size
 
; !p.background=65535

  xwin = 720
  ywin = 720

; scale the window for the postscript plot

  xps  = 6.
  yps  = xps * ywin / xwin

; define the universal time and longitude

  hr   = time(1,ntm)
  phr  = string(format='(i2)',hr)
  while(((i=strpos(phr,' '))) ne -1) do strput,phr,'0',i
  min  = time(2,ntm)
  pmin = string(format='(i2)',min)
  while(((i=strpos(pmin,' '))) ne -1) do strput,pmin,'0',i
  t    = 'UT: '+phr+pmin

  day  = nday + fix(time(4,ntm)/24.0)
  pday = 'Day ' + string(format='(i3)',day)
  
  lon00=0.
  lon0 =  lon00
  
; define the local time and longitude

  sec  = time(3,ntm)
  utsec  = hr * 3600. + min * 60. + sec
  uthr   = utsec / 3600.
  lthr   = uthr + lon0/15.
  if lthr ge 24. then lthr = lthr - 24.

;  lon00  = (((uthr)*15.)mod 360.+lon0)mod 360.

;  lon00  = -(uthr*15.)mod 360.
; lon00=0.
  print,uthr,lon00

; set view

  lat0 =  90.
  lon0 =  lon00
  lmt  = [0,0,90,360]


;  lon1  = fltarr(nf,nl+1)
;  lat1  = fltarr(nf,nl+1)
;  fcp1  = fltarr(nf,nl+1)

;  lon1  = reform(blon(nz-1,0:nf-1,*))
;  lat1(0:nf-2,*)  = reform(blat(nz-1,0:nf-2,*))
;  lat1(nf-1,*)    = 90.
;  fcp1(0:nf-2,*)  = reform(transpose(phi(*,0:nf-2,ntm)))*300./1.e3
;  fcp1(nf-1,*)    = 0.
;  phi0 = (max(fcp1)-min(fcp1))/2.
;  fcp1 = fcp1 + phi0

  lon1  = reform(blon(nz-1,0:nf-2,*))
  lat1  = reform(blat(nz-1,0:nf-2,*))
  fcp1  = reform(transpose(phi(*,0:nf-2,ntm)))*300./1.e3

;  fcp1  = reform(phip(0:nf-2,*,ntm))*300./1.e3

;;  phi0 = (max(fcp1)-min(fcp1))/2.
;;  fcp1 = fcp1 + phi0

;  fcp1  = fcp1 - max(fcp1) ; max(fcp1) > 0

  colorphi = 0

  window,xsize=xwin,ysize=ywin,iwin

 latmax  = max(lat1)
 latmin  = min(lat1)
 lonmax  = max(lon1)
 lonmin  = min(lon1)
 fmax1    = max(fcp1)*1.05
 fmin1    = min(fcp1)*.9

 print,fmin1,fmax1
 print,latmin,latmax
 print,lonmin,lonmax

 fmin1   = -100.
 fmax1   = 100.

; fmin1   = -phi0 
; fmax1   =  phi0 

; fmin1  = -60.
; fmax1  =  60.

 print,fmin1,fmax1


 nsteps  = 40
 step    = ( fmax1 - fmin1 ) / (nsteps-1.)
 levsf1  = findgen(nsteps) * step + fmin1

 step    = ( fmax1 - fmin1 ) / (nsteps-1.)
 levsf2  = findgen(nsteps) * step + fmin1
 levsf2  = levsf1


 !p.thick=2
 map_set,/orthographic,/grid,/label,/isotropic,lat0,lon0,$
   limit=lmt,/noerase
 contour,fcp1,lon1,lat1,levels=levsf1,/overplot,/cell_fill
 contour,fcp1,lon1,lat1,levels=levsf2,/overplot,color=colorphi
 map_set,/orthographic,/grid,/label,/isotropic,lat0,lon0,$
         /noerase,color=255;,limit=lmt

; colorbar

  posc = [.65,.05,.95,.05] 
  xszw = ( posc(2) - posc(0) ) * !D.X_Vsize
  bar  = bindgen(256) # replicate(1b,20)
  tvscl,congrid(bar,xszw,20,/interp),posc(0),posc(1),/normal   
  cxs = fltarr(3)
  cxs = [fmin1,.5*(fmin1+fmax1),fmax1]
  cxsl= [string(format='(f5.1)',cxs(0)),$
         string(format='(f5.1)',cxs(1)),$
         string(format='(f5.1)',cxs(2))]
  cxsl= [string(format='(e8.1)',cxs(0)),$
         string(format='(e8.1)',cxs(1)),$
         string(format='(e8.1)',cxs(2))]
  plot,cxs,/noerase,/nodata,position=posc,xticks=2,$
    ystyle=4,xstyle=1,xtickname=cxsl,charsize=1.2 ,col=255


  lhr    = ceil(lthr) - 1
  print,'   hUT = ',hrut(ntm)
  plhr  = string(format='(i2)',lhr)
  while(((i=strpos(plhr,' '))) ne -1) do strput,plhr,'0',i
  lmin  = ceil ( ( lthr - lhr ) * 60. ) - 1
  plmin = string(format='(i2)',lmin)
  while(((i=strpos(plmin,' '))) ne -1) do strput,plmin,'0',i
  plt    = 'LT: '+plhr+plmin
  

; print the universal time and local time

  xp = .1
  yp = .9
  xyouts,xp,yp,pday,size=1.4,color=255,/normal
  xp = .1
  yp = .87
  xyouts,xp,yp,t,size=1.4,color=255,/normal


;stop

; line plot of density along field line

;  !P.Position=[.65,.75,.85,.925]
;  jlon = min(where ( glon(0,*,nf) ge lon0-2 and glon(0,*,nf) le lon0+2,count ))
;  plot,value(*,jlon),zalt(*,jlon,nf),/noerase, $
;    xrange=[valmin,valmax],yrange=[100,500], $
;    xtitle = '!6Electron Density (log !8n!de!n!6)', $
;    ytitle = '!6Altitude along B (km)'
;  !P.Position=0
;stop
; -----------------------
; set plot to postscript
; -----------------------

  postscript = 1
 
  if postscript eq 1 then begin

  !p.charthick = 4
  !p.thick     = 4
  xyth         = 4

  set_plot,'ps'
  device,file='fig.ps',bits_per_pixel=8,/color, $
         xsize=xps,ysize=yps,/inches,xoffset=1.25


 map_set,/orthographic,/grid,/label,/isotropic,lat0,lon0,$
   limit=lmt,/noerase
 contour,fcp1,lon1,lat1,levels=levsf1,/overplot,/cell_fill
 contour,fcp1,lon1,lat1,levels=levsf2,/overplot,color=colorphi
 map_set,/orthographic,/grid,/label,/isotropic,lat0,lon0,$
         /noerase,color=255;,limit=lmt

; colorbar

  posc = [.62,.95,.92,.95]   

  xszc = ( posc(2) - posc(0) ) * xps
  yszc = ( posc(3) - posc(1) ) * yps
  xstc = posc(0) * xps
  ystc = posc(1) * yps
  xszw = ( posc(2) - posc(0) ) * xwin

  bar  = bindgen(256) # replicate(1b,10)
  tvscl,congrid(bar,xwin,10,/interp),$
        xstc,ystc,xsize=xszc,ysize=.2,/inches
  cxs = fltarr(3)
  cxs = [fmin1,.5*(fmin1+fmax1),fmax1]
  cxsl= [string(format='(f5.1)',cxs(0)),$
         string(format='(f5.1)',cxs(1)),$
         string(format='(f5.1)',cxs(2))]
  cxsl= [string(format='(e8.1)',cxs(0)),$
         string(format='(e8.1)',cxs(1)),$
         string(format='(e8.1)',cxs(2))]
  plot,cxs,/noerase,/nodata,position=posc,xticks=2,$ 
       ystyle=4,xstyle=1,charsize=1.,$
       xtickname=cxsl,color=0,$
       xtitle=cbt

  xp = .1
  yp = .95
  xyouts,xp,yp,pday,size=1.2,color=0,/normal
  xp = .1
  yp = .9
  xyouts,xp,yp,t,size=1.2,color=0,/normal

  device,/close

  endif

; reset the display to the screen

  set_plot,'x'

  !p.charthick = 1
  !p.thick     = 1

 end

