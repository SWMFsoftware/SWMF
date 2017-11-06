;-------------------------------------------------------------------------------
;
;                                plot_drift.pro
; 
; IDL procedure reads *.fls file and plot the pitch-angle averaged flux,
; pressure, density and pitch-angle distribution at the equatorial plane.
;
; Created on 18 September 2007 by Mei-Ching Fok, Code 673, NASA GSFC.
; 
; Modification History
; May 28, 2015
;   * Edited to plot the calculated drifts of (relatavistic) particles.
; June 8, 2008
;   * Add the calculation of temperature
; June 10, 2008
;   * The pressure calculation is valid in relativitic limit.
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
pro plot_equator,var,vmax,vmin,xi,bar_lab,halfl,ro,mlto,lat,x_wsize,y_wsize, $
                 dir,blackc,g_thick,iround
;-------------------------------------------------------------------------------

; Make levels
nlevel=59
lvl=fltarr(nlevel)    &    colr=intarr(nlevel)
dlvl=(vmax-vmin)/(nlevel-1)
colrmax=254          
colrmin=1
ncolor=colrmax-colrmin+1
dcolr=(float(colrmax)-float(colrmin))/(nlevel-1)
for i=0,nlevel-1 do begin
    lvl(i)=vmin+i*dlvl
    colr(i)=round(float(colrmin)+i*dcolr)
endfor
yi=0.30
plt_size=y_wsize*0.5
xf=xi+plt_size/x_wsize
yf=yi+plt_size/y_wsize

; force var equal or gt vmin 
if (min(var) lt vmin) then var(where(var le vmin))=vmin

; polyfill background color in black
polyfill,[xi,xf,xf,xi,xi],[yi,yi,yf,yf,yi],color=blackc,/normal

; plot flux or PA anisotropy
if (dir eq 1) then phi=mlto               ; Sun to the left
if (dir eq 2) then phi=mlto+!pi           ; Sun to the right
polar_contour,var,phi,ro,xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
      levels=lvl,c_colors=colr,xstyle=1,ystyle=1,xticks=4,yticks=4, $
      xtitle='RE',pos=[xi,yi,xf,yf],color=blackc,/fill,/noerase
dim=size(ro)
for i=0,dim(1)-1 do begin
    if (lat(i) ge 58.) then begin
       ihigh=i
;      print,'ihigh = ',ihigh
       goto,continue
    endif
endfor
continue:
;for j=0,dim(2)-1 do oplot,ro(0:dim(1)-1,j),phi(0:dim(1)-1,j),/polar
;for i=0,dim(1)-1 do oplot,ro(i,0:dim(2)-1),phi(i,0:dim(2)-1),/polar
;oplot,ro(ihigh,0:dim(2)-1),phi(ihigh,0:dim(2)-1),/polar,color=254

; Add color bar and label
xb=xi
yb=0.090
dx=0.006
y1=yb+0.05
for i=0,nlevel-1 do begin
    x1=xb+i*dx
    x2=x1+dx*1.02
    polyfill,[x1,x2,x2,x1,x1],[yb,yb,y1,y1,yb],color=colr(i),/normal
endfor
v_min=string(vmin,'(f6.1)')
if (vmax lt 100.) then v_max=string(vmax,'(f6.1)')
if (vmax ge 100.) then v_max=string(vmax,'(e9.1)')
;if (bar_lab eq 'pitch angle anisotropy') then v_min='field aligned'
;if (bar_lab eq 'pitch angle anisotropy') then v_max='perpendicular'
if (bar_lab eq 'pitch angle anisotropy') then v_min='-1.0' 
if (bar_lab eq 'pitch angle anisotropy') then v_max='1.0'
xyouts,xb,0.6*yb,v_min,alignment=0.5,size=1.2,color=blackc,/normal
xyouts,x2,0.6*yb,v_max,alignment=0.5,size=1.2,color=blackc,/normal
if (bar_lab eq 'pitch angle anisotropy') then begin
   xyouts,xb,0.20*yb,'field aligned',alignment=0.5,size=1.2,color=blackc,/normal
   xyouts,x2,0.20*yb,'perpendicular',alignment=0.5,size=1.2,color=blackc,/normal
endif
xyouts,0.5*(x2+xb),1.7*yb,bar_lab,alignment=0.5,size=1.1,color=blackc,/normal

; fill cirle of radius ro(0,0) centered the earth with color2
color2=colrmin
if (bar_lab eq 'pitch angle anisotropy') then begin
   m=(nlevel-1)/2
   color2=colr(m)
endif

; Draw earth , geosynchronous, round background if iround=1
npt=400                             ; no. of points in a circle
npt2=npt/2
night_x=fltarr(npt2+2)     &    day_x=night_x
night_y=fltarr(npt2+2)     &    day_y=night_y
e1x=fltarr(npt+1)    &   e1y=e1x
del = 2.0 * !pi / npt
for i = 0,npt do begin
    i1 = i - npt2
    if (dir eq 1) then ang = float(i) * del + !pi/2. 
    if (dir eq 2) then ang = float(i) * del - !pi/2.
    cosa = cos(ang)
    sina = sin(ang)
    if (i le npt2) then day_x(i) = cosa
    if (i le npt2) then day_y(i) = sina
    if (i ge npt2) then night_x(i1) = cosa
    if (i ge npt2) then night_y(i1) = sina
    e1x(i) = 6.6 * cosa
    e1y(i) = 6.6 * sina
endfor
day_x(npt2+1) = day_x(0)
day_y(npt2+1) = day_y(0)
night_x(npt2+1) = night_x(0)
night_y(npt2+1) = night_y(0)
polyfill,ro(0,0)*night_x,ro(0,0)*night_y,color=color2
polyfill,ro(0,0)*day_x,ro(0,0)*day_y,color=color2
polyfill,night_x,night_y,color=blackc
polyfill,day_x,day_y,color=255
if (iround eq 1) then begin
   e2x=fltarr(npt+6)   &   e2y=e2x
   e2x(0:npt)=e1x(0:npt)*halfl/6.6
   e2y(0:npt)=e1y(0:npt)*halfl/6.6
   e2x(npt+1:npt+5)=[-halfl,-halfl,halfl,halfl,0]
   if (dir eq 1) then e2y(npt+1:npt+5)=[halfl,-halfl,-halfl,halfl,halfl]
   if (dir eq 2) then e2y(npt+1:npt+5)=[-halfl,halfl,halfl,-halfl,-halfl]
   polyfill,e2x,e2y,color=blackc
endif
oplot,e1x,e1y,color=255,thick=g_thick     ; plot geosynchronous

return
end


;-----------------------------------------------------------------------------
; main routine
;-----------------------------------------------------------------------------

close,/all
test_string=''
; get colors for color table palette
read,'color table, 1=rainbow with white, 2=HENA-type => ',icr
if (icr eq 1) then begin
   loadct,39       ;  rainbow with white
   tvlct, red, green, blue, /get
   blackc=0
endif
if (icr ne 1) then begin           ; color table of HENA plots
   red=fltarr(256)    &    green=fltarr(256)    &    blue=fltarr(256)
   openr,1,'hena.tbl'
   readf,1,red,green,blue
   close,1
   tvlct, red, green, blue
   blackc=254
endif
!p.color=blackc
!p.charsize=1.5

; read file name and energy, sina, lat information 
fhead=' '
read,'enter the file head (e.g., 2000_225_h) => ',fhead
openr,2,fhead+'.vl'
openr,3,fhead+'.vp'
readf,2,rc,ir,ip,je,ig   ; je=number of energy grid, ig=number of pa grid
readf,3,rc,ir,ip,je,ig   ; je=number of energy grid, ig=number of pa grid
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
energy=fltarr(je)    &   sina=fltarr(ig)   &   lat=fltarr(ir)   &  ece=energy
Ebound=fltarr(je+1)
readf,2,energy
readf,3,energy
readf,2,sina
readf,3,sina
PA_bins_deg=asin(sina)*!radeg
readf,2,lat
readf,3,lat
read,'number of data set in time => ',ntime
read,'scale?  0 = linear scale, or 1 = log scale => ',ilog
read,'Sun to the?   1 = left,  2 = right => ',dir
read,'half length of the plot => ', halfl
read,'shape of the plot?  0 = square,  1 =  round => ',iround

; Calculate Ebound    
;; for k=1,je-1 do Ebound(k)=sqrt(energy(k-1)*energy(k))
;; Ebound(0)=energy(0)*energy(0)/Ebound(1)
;; Ebound(je)=energy(je-1)*energy(je-1)/Ebound(je-1)
;; if (Ebound(je) lt 999.5) then ilab='(i3.3)'
;; if (Ebound(je) ge 999.5) then ilab='(i4.4)'
;; if (Ebound(je) lt 999.5) then flab='(f5.1)'
;; if (Ebound(je) ge 999.5) then flab='(f6.1)'
 
; setup energy and species
nsp=strlen(fhead)-2
spe=strmid(fhead,nsp,2)                ; species
if (spe eq '_h') then species='H+'
if (spe eq '_s') then species='SW H+'
if (spe eq 'Po') then species='PolarWind H+'
if (spe eq 'Pl') then species='PlasmasphericWind H+'
if (spe eq '_o') then species='O+'
IF (Spe eq '_e') then species='e-'

; setup arrays
vl=fltarr(ntime,nr,nmlt+1,je,ig)
vp=fltarr(ntime,nr,nmlt+1,je,ig)
roa=fltarr(ntime,nr,nmlt+1)   &  mltoa=roa   
irm=intarr(ntime,nmlt)        &  iba=irm
ro=fltarr(nr,nmlt+1)          &  mlto=ro
vlplot=fltarr(nr,nmlt+1)
vpplot=fltarr(nr,nmlt+1)         
houra=fltarr(ntime)
vl_temp=fltarr(ig)
vp_temp=fltarr(ig)
;; cosa(*)=cos(asin(sina(*)))
;; for m=0,ig-1 do begin
;;     if (m eq 0) then sina0=0.
;;     if (m gt 0) then sina0=0.5*(sina(m)+sina(m-1))
;;     if (m eq (ig-1)) then sina1=1.
;;     if (m lt (ig-1)) then sina1=0.5*(sina(m)+sina(m+1))
;;     dmu(m)=sqrt(1.-sina0*sina0)-sqrt(1.-sina1*sina1)
;; endfor
;; dE=fltarr(je)
;; for k=0,je-1 do dE(k)=Ebound(k+1)-Ebound(k)

; factors for calculating drift velocity
velf=1                          ; velocity from rad/s to km/s

; Calculate ece, (E+Eo)/c/sqrt(E(E+2Eo)). c in cm/s
if (spe eq '_h' or spe eq '_s' or spe eq 'Po' or spe eq 'Pl') then mass=1.
if (spe eq '_o') then mass=16.
if (spe eq '_e') then mass=5.4462e-4
Eo=mass*1.673e-27*3.e8*3.e8/1.6e-16        ; rest-mass energy in keV
;for k=0,je-1 do ece(k)=(energy(k)+Eo)/3.e10/sqrt(energy(k)*(energy(k)+2.*Eo))

; Read drift velocities
for n=0,ntime-1 do begin
    readf,2,hour
    readf,3,hour
    houra(n)=hour
    print,n,hour
    for i=0,nr-1 do begin
        for j=0,nmlt-1 do begin
           readf,2,lat1,mlt,ro1,mlto1,bo ;,irm1,iba1
           readf,3,lat1,mlt,ro1,mlto1,bo ;,irm1,iba1
           ;; irm(n,j)=irm1
           ;; iba(n,j)=iba1
           roa(n,i,j)=ro1
           mltoa(n,i,j)=mlto1
           for k=0,je-1 do begin
              ;; for m=0,ig-1 do begin
                 readf,2,vl_temp
                 readf,3,vp_temp

                 T=energy(k)/Eo
                 gamma=T+1

                 ;; if $
                 ;;    ((ro1 ge 6.) and $
                 ;;     (ro1 le 7.) and $
                 ;;     (k eq je-1)) $
                 ;; then begin
                 ;;    print,'lat1,mlt,ro1,mlto1,bo: ', $
                 ;;          lat1,mlt,ro1,mlto1,bo
                 ;;    print,'n,i,j,k,m: ',n,i,j,k,m
                 ;;    print,'vp_temp: ',vp_temp(11)
                 ;;    print,'ro1*vp_temp*velf', $
                 ;;          ro1*vp_temp(11)*velf
                 ;;    print,'gamma*ro1*vp_temp*velf', $
                 ;;          gamma*ro1*vp_temp(11)*velf
                 ;;    ;; read,test_string
                 ;; endif

                 vl(n,i,j,k,*)=ro1*vl_temp*velf
                 vp(n,i,j,k,*)=vp_temp
;                 vp(n,i,j,k,*)=ro1*vp_temp*velf
              
           endfor
           
        endfor
        
     endfor
            
    mltoa(n,*,*)=mltoa(n,*,*)*!pi/12. ; mlto in radian
    ; periodic boundary condition
    roa(n,*,nmlt)=roa(n,*,0)
    mltoa(n,*,nmlt)=mltoa(n,*,0)
    vl(n,*,nmlt,*,*)=vl(n,*,0,*,*)
    vp(n,*,nmlt,*,*)=vp(n,*,0,*,*)
 endfor
close,2

new_plot:

!p.charthick=1.5
g_thick=1.5

plot_continue='y'

while (plot_continue eq 'y') do begin
; choose which energy to be displaced
print,'energy bins (keV): '
for i=0,je-1 do print,i,Energy(i),format='("  (",i2,") ",f7.2)'
read, 'Which energy bin (ie)? => ',ie
e0=energy(ie)
;; e1=Ebound(ie+1)
print,'PA bins (Deg.): '
for i=0,ig-1 do print,i,PA_bins_deg(i),format='("  (",i2,") ",f4.1)'
read, 'Which PA bin (iPA)? => ',iPA
PA_sel=PA_bins_deg(iPA)

yon=' '
; Determines max and min values for the radial velocity
vlmax=max(vl(*,*,*,ie,iPA))    
print,' vlmax = ',vlmax
read,'Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,'enter new vlmax => ',vlmax
vlmin=min(vl(*,*,*,ie,iPA))
if (ilog eq 1) then vlmin=vlmax-3.   
print,' vlmin = ',vlmin
read,'Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,'enter new vlmin => ',vlmin

;; endelse

; Determines max and min values for the poloidal velocity
vpmax=max(vp(*,*,*,ie,iPA))    
print,' vpmax = ',vpmax
read,'Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,'enter new vpmax => ',vpmax
vpmin=min(vp(*,*,*,ie,iPA))
if (ilog eq 1) then vpmin=vpmax-3.   
print,' vpmin = ',vpmin
read,'Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,'enter new vpmin => ',vpmin

; smooth data
read,'Do you want to smooth the data? (y/n) => ',yon
nframe=0
if (yon eq 'y') then read,'enter number of intermediate frame => ',nframe
ntnew=(ntime-1)*(nframe+1)+1
ronew=fltarr(ntnew,nr,nmlt+1)  &  mltonew=ronew  
vlnew=ronew
vpnew=ronew
if (nframe eq 0) then begin
   hournew=houra        
   ronew=roa        
   mltonew=mltoa   

   vlnew=vl(*,*,*,ie,iPA)
   vpnew=vp(*,*,*,ie,iPA)
   
endif else begin
   hournew=interpol(houra,ntnew)
   for i=0,nr-1 do begin 
       for j=0,nmlt do begin
           ronew(*,i,j)=interpol(roa(*,i,j),ntnew)
           mltonew(*,i,j)=interpol(mltoa(*,i,j),ntnew)
           vlnew(*,i,j)=interpol(vl(*,i,j,ie,iPA),ntnew)
           vpnew(*,i,j)=interpol(vp(*,i,j,ie,iPA),ntnew)
       endfor
   endfor
endelse

vptest_check=''
read,'Plot vptest? (y/n) => ',vptest_check

; Setup window
;; ips=0
;; set_plot,'x'
x_wsize=684
y_wsize=480
;; window,1,xpos=200,ypos=100,xsize=x_wsize,ysize=y_wsize

; Plot radial and poloidal velocities
for n=0,ntnew-1 do begin
   setdevice, $
      'Drifts'+spe+'_E='+STRTRIM(STRING(energy(ie),FORMAT='(F07.2)'),2)+ $
      ',PA='+STRTRIM(STRING(PA_sel,FORMAT='(I02)'),2)+'.eps',/eps
   polyfill,[0.,1.,1.,0.,0.],[0.,0.,1.,1.,0.],color=255,/normal ; white BG
    hour=hournew(n)
    ro(*,*)=ronew(n,*,*)
    mlto(*,*)=mltonew(n,*,*)
    vlplot(*,*)=vlnew(n,*,*)
    vpplot(*,*)=vpnew(n,*,*)

; plot poloidal drift
    xi=0.56
    bar_lab='Poloidal Drift (km/s)'
    if (ilog eq 1) then bar_lab='log '+bar_lab
    plot_equator,vpplot,vpmax,vpmin,xi,bar_lab,halfl,ro,mlto,lat, $
                 x_wsize,y_wsize,dir,blackc,g_thick,iround

; plot radial drift
    xi=0.07
    bar_lab='Radial Drift (km/s)'
    ;; if (vptest_check eq 'y') then bar_lab='Poloidal Drift Test (km/s)'
    if (ilog eq 1) then bar_lab='log '+bar_lab
    plot_equator,vlplot,vlmax,vlmin,xi,bar_lab,halfl,ro,mlto,lat, $
                 x_wsize,y_wsize,dir,blackc,g_thick,iround

;label the plot
    hours=string(fix(hour),'(i3)')
    minute=round((hour-float(fix(hour)))*60.)
    minutes=string(minute,'(i2.2)')
    elabel=string(e0);+' - '+string(e1,flab)
    
    xyouts,0.07,0.93,fhead,size=2.,color=blackc,/normal
    xyouts,0.5,0.93,hours+':'+minutes+' UT',size=2.,alignment=0.5, $
           color=blackc,/normal 
    xyouts,0.5,0.85,elabel+' keV, PA='+ $
           strtrim(string(PA_sel,FORMAT="(F4.1)"),2)+' deg. '+species, $
           size=2.,alignment=0.5,color=blackc,/normal
    closedevice

    if (vptest_check eq 'y') then begin
       
       vptest=fltarr(nr,nmlt+1)

       T=energy(ie)/Eo
       beta=sqrt(T*(T+2)/(T+1)^2)
       gamma=T+1
       W=(energy(ie)+Eo)*1.602e-16	; keV to J
       q=1.602e-19                      ; Coulomb
       B_E=30574e-9                     ; Tesla
       Bo=B_E                           ; Tesla
       R_E=6.378e6                      ; meters
       m_i=1.67e-27                     ; kg

       ;; if $
       ;;    (spe eq '_h' or spe eq '_s' or spe eq 'Po' or spe eq 'Pl') $
       ;; then $
       ;;    c_d=8.481
       ;; if (spe eq '_o') then c_d=8.481/16.
       ;; if (spe eq '_e') then c_d=1.557e4
       c_d=2.*!pi*q*B_E*R_E^2/mass/m_i/3e8/3e8
       
       t_drift_geo=c_d/6.6/gamma/beta/beta*(1-1/3.*sina(ipa)^0.62)
       vdrift_geo=2*!pi*6.6*R_E/t_drift_geo
       if (spe ne '_e') then vdrift_geo=-vdrift_geo
       print,'GEO Azimuthal drift time (s): ',t_drift_geo
       print,'GEO Azimuthal drift time (km/s): ',vdrift_geo/1e3
       
       ieqn=0
       ;; read,'Drift Equation?  0 = Baumjohann, or 1 = Walt => ',ieqn
       for iL=0,nr-1 do begin
             
          for jMLT=0,nmlt-1 do begin

;;              if ieqn eq 0 then $
                
;; ; Equation from Baumjohann; Does not explicitly account for
;; ; relativistic effects
;;                 vptest(iL,jMLT)= $
;;                 6.*ro(iL,jMLT)^2*W/q/B_E/R_E* $
;;                 (0.35+0.15*sina(iPA))*1e-3 $ ; in km/s
;;              else $
; Equation from Walt
             t_drift= $
                c_d/ro(iL,jMLT)/gamma/beta/beta*(1-1/3.*sina(ipa)^0.62)
             vptest(iL,jMLT)= $
                2*!pi*ro(iL,jMLT)*R_E/t_drift/1e3
             
          endfor
          
       endfor

       vptest(*,nMLT)=vptest(*,0)
       
       if $
          (spe ne '_e') $
       then begin
          
          vptest=-vptest

          vptestmin=vdrift_geo/1e3 ;min([max(vptest),max(vpplot)])
          vptestmax=0.

       endif else begin

; Determines max and min values for the test poloidal velocity
          vptestmax=vdrift_geo/1e3 ;min([max(vptest),max(vpplot)])
          vptestmin=min([min(vptest),min(vpplot)])

       endelse

       setdevice, $
          'Drift_Comparison'+spe+'_E='+ $
          STRTRIM(STRING(energy(ie),FORMAT='(F07.2)'),2)+',PA='+ $
          STRTRIM(STRING(PA_sel,FORMAT='(I02)'),2)+'.eps',/eps
       polyfill,[0.,1.,1.,0.,0.],[0.,0.,1.,1.,0.],color=255,/normal ; white BG

; plot poloidal drift
       xi=0.56
       bar_lab='Poloidal Drift (km/s)'
       if (ilog eq 1) then bar_lab='log '+bar_lab
       bar_lab='Model '+bar_lab
       plot_equator,vpplot,vptestmax,vptestmin,xi,bar_lab,halfl,ro,mlto,lat, $
                    x_wsize,y_wsize,dir,blackc,g_thick,iround
       
; plot radial drift
       xi=0.07
       bar_lab='Poloidal Drift (km/s)'
       if (ilog eq 1) then bar_lab='log '+bar_lab
       ;; if (ieqn eq 0) then bar_lab='Baumjohann '+bar_lab
       ;; if (ieqn eq 1) then bar_lab='Walt '+bar_lab
       bar_lab='Analytical '+bar_lab
       plot_equator,vptest,vptestmax,vptestmin,xi,bar_lab,halfl,ro,mlto,lat, $
                    x_wsize,y_wsize,dir,blackc,g_thick,iround

;label the plot
       hours=string(fix(hour),'(i3)')
       minute=round((hour-float(fix(hour)))*60.)
       minutes=string(minute,'(i2.2)')
       elabel=string(e0);+' - '+string(e1,flab)
       
       xyouts,0.07,0.93,fhead,size=2.,color=blackc,/normal
       xyouts,0.5,0.93,hours+':'+minutes+' UT',size=2.,alignment=0.5, $
              color=blackc,/normal 
       xyouts,0.5,0.85,elabel+' keV, PA='+ $
              strtrim(string(PA_sel,FORMAT="(F4.1)"),2)+' deg. '+species, $
              size=2.,alignment=0.5,color=blackc,/normal
       closedevice
       
       setdevice, $
          'Drift_Percent_Diff'+spe+'_E='+ $
          STRTRIM(STRING(energy(ie),FORMAT='(F07.2)'),2)+',PA='+ $
          STRTRIM(STRING(PA_sel,FORMAT='(I02)'),2)+'.eps',/eps

       polyfill,[0.,1.,1.,0.,0.],[0.,0.,1.,1.,0.],color=255,/normal ; white BG

; plot poloidal drift
       xi=0.56
       bar_lab='Absolute Percent Difference'
       plot_equator,abs((vpplot-vptest)/vptest)*100., $
                    100.,0.,xi,bar_lab,halfl,ro,mlto,lat, $
                    x_wsize,y_wsize,dir,blackc,g_thick,iround

;label the plot
       hours=string(fix(hour),'(i3)')
       minute=round((hour-float(fix(hour)))*60.)
       minutes=string(minute,'(i2.2)')
       elabel=string(e0);+' - '+string(e1,flab)
       
       xyouts,0.07,0.93,fhead,size=2.,color=blackc,/normal
       xyouts,0.5,0.93,hours+':'+minutes+' UT',size=2.,alignment=0.5, $
              color=blackc,/normal 
       xyouts,0.5,0.85,elabel+' keV, PA='+ $
              strtrim(string(PA_sel,FORMAT="(F4.1)"),2)+' deg. '+species, $
              size=2.,alignment=0.5,color=blackc,/normal
       closedevice
       
    endif
    
 endfor

read,'Do you want to continue making plots? (y/n) => ',plot_continue
if plot_continue eq 'y' then print,'Previous energy index: ',ie
if plot_continue eq 'y' then print,'Previous PA index: ',iPA

endwhile
    
end
