;-------------------------------------------------------------------------------
;
;                                plot_fls.pro
; 
; IDL procedure reads *.fls file and plot the pitch-angle averaged flux,
; pressure, density and pitch-angle distribution at the equatorial plane.
;
; Created on 18 September 2007 by Mei-Ching Fok, Code 673, NASA GSFC.
; 
; Modification History
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
v_min=string(vmin,'(f4.1)')
if (vmax lt 100.) then v_max=string(vmax,'(f4.1)')
if (vmax ge 100.) then v_max=string(vmax,'(e7.1)')
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
openr,2,fhead+'.fls'
readf,2,rc,ir,ip,je,ig   ; je=number of energy grid, ig=number of pa grid
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
energy=fltarr(je)    &   sina=fltarr(ig)   &   lat=fltarr(ir)   &  ece=energy
Ebound=fltarr(je+1)
readf,2,energy
readf,2,sina
readf,2,lat
read,'number of data set in time => ',ntime
read,'scale?  0 = linear scale, or 1 = log scale => ',ilog
read,'Sun to the?   1 = left,  2 = right => ',dir
read,'half length of the plot => ', halfl
read,'shape of the plot?  0 = square,  1 =  round => ',iround

; Calculate Ebound    
for k=1,je-1 do Ebound(k)=sqrt(energy(k-1)*energy(k))
Ebound(0)=energy(0)*energy(0)/Ebound(1)
Ebound(je)=energy(je-1)*energy(je-1)/Ebound(je-1)
if (Ebound(je) lt 999.5) then ilab='(i3.3)'
if (Ebound(je) ge 999.5) then ilab='(i4.4)'
if (Ebound(je) lt 999.5) then flab='(f5.1)'
if (Ebound(je) ge 999.5) then flab='(f6.1)'
 
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
plsfluxk=fltarr(ntime,nr,nmlt+1,je)    &  anisok=plsfluxk
plsfluxa=fltarr(ntime,nr,nmlt+1)       &  anisoa=plsfluxa
roa=fltarr(ntime,nr,nmlt+1)   &  mltoa=roa   
irm=intarr(ntime,nmlt)        &  iba=irm
ro=fltarr(nr,nmlt+1)          &  mlto=ro     &   plsflux=ro      &  aniso=ro
houra=fltarr(ntime)
fy=fltarr(ig)      &      ypar=fy      &     dmu=fy        &    cosa=fy
cosa(*)=cos(asin(sina(*)))
for m=0,ig-1 do begin
    if (m eq 0) then sina0=0.
    if (m gt 0) then sina0=0.5*(sina(m)+sina(m-1))
    if (m eq (ig-1)) then sina1=1.
    if (m lt (ig-1)) then sina1=0.5*(sina(m)+sina(m+1))
    dmu(m)=sqrt(1.-sina0*sina0)-sqrt(1.-sina1*sina1)
endfor
dE=fltarr(je)
for k=0,je-1 do dE(k)=Ebound(k+1)-Ebound(k)

; factor of calculating ion pressure, density and temperature
pf=8.*!pi/3.*1.e6*1.6e-16*1.e9    ; pressure in nPa
nf=4.*!pi                         ; density in cm-3
Tf=2./3.                          ; T = 2/3 * mean E

; Calculate ece, (E+Eo)/c/sqrt(E(E+2Eo)). c in cm/s
if (spe eq '_h' or spe eq '_s' or spe eq 'Po' or spe eq 'Pl') then mass=1.
if (spe eq '_o') then mass=16.
if (spe eq '_e') then mass=5.4462e-4
Eo=mass*1.673e-27*3.e8*3.e8/1.6e-16        ; rest-mass energy in keV
for k=0,je-1 do ece(k)=(energy(k)+Eo)/3.e10/sqrt(energy(k)*(energy(k)+2.*Eo))

; Read fluxes and calculate PA anisotropy
for n=0,ntime-1 do begin
    readf,2,hour
    houra(n)=hour
    print,n,hour
    for i=0,nr-1 do begin
        for j=0,nmlt-1 do begin
            readf,2,lat1,mlt,ro1,mlto1,bo  ;,irm1,iba1
;           irm(n,j)=irm1
          ; iba(n,j)=iba1
            roa(n,i,j)=ro1
            mltoa(n,i,j)=mlto1
            for k=0,je-1 do begin
                readf,2,fy
                ; calculate pitch-angle (PA) averaged flux and PA anisotropy
                plsfluxk(n,i,j,k)=0.
                anisok(n,i,j,k)=0.
                fpara=0.
                fperp=0. 
                for m=0,ig-1 do begin
                    fydmu=fy(m)*dmu(m)
                    plsfluxk(n,i,j,k)=plsfluxk(n,i,j,k)+fydmu
                    fpara=fpara+fydmu*cosa(m)*cosa(m)
                    fperp=fperp+fydmu*sina(m)*sina(m)/2.
                endfor
                ft=fpara+fperp
                if (ft gt 0.) then anisok(n,i,j,k)=(fperp-fpara)/ft
            endfor
        endfor
    endfor
    mltoa(n,*,*)=mltoa(n,*,*)*!pi/12.                     ; mlto in radian
    ; periodic boundary condition
    roa(n,*,nmlt)=roa(n,*,0)
    mltoa(n,*,nmlt)=mltoa(n,*,0)
    plsfluxk(n,*,nmlt,*)=plsfluxk(n,*,0,*)
    anisok(n,*,nmlt,*)=anisok(n,*,0,*)
endfor
close,2

new_plot:

!p.charthick=1.5
g_thick=1.5

; choose which energy to be displaced
print,'energy bin (keV): '
for i=0,je-1 do print,i,Ebound(i:i+1),format='("  (",i2,") ",f7.2," - ",f7.2)'
read, 'lower and upper energy bins (ie1,ie2) => ',ie1,ie2
read,'plot? (1)flux,  (2)pressure,  (3)density,  (4)temperature => ',iopt
iunit=0
if (iopt eq 1) then $
    read,'unit? (1)flux/(keV cm2 s sr), (2)flux/(cm2 s sr) => ',iunit
e0=Ebound(ie1)
e1=Ebound(ie2+1)
fmiddle=string(round(e0),ilab)+'-'+string(round(e1),ilab)+'keV'
if (iopt eq 1) then fmiddle=fmiddle+'_flux' 
if (iopt eq 2) then fmiddle=fmiddle+'_pressure'
if (iopt eq 3) then fmiddle=fmiddle+'_density'
if (iopt eq 4) then fmiddle=fmiddle+'_temperature'

; calculate energy-integrated flux or pressure or density or temperature
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        for j=0,nmlt do begin
            plsfluxa(n,i,j)=0.
            anisoa(n,i,j)=0.
            weiSum=0.
            denSum=0.
            for k=ie1,ie2 do begin    
                if (iopt eq 1) then wei=dE(k)
                if (iopt eq 2 or iopt eq 4) then wei=dE(k)*energy(k)*ece(k)
                if (iopt eq 3) then wei=dE(k)*ece(k)           
                weiSum=weiSum+wei
                if (iopt eq 4) then denSum=denSum+plsfluxk(n,i,j,k)*dE(k)*ece(k)
                plsfluxa(n,i,j)=plsfluxa(n,i,j)+plsfluxk(n,i,j,k)*wei  
                anisoa(n,i,j)=anisoa(n,i,j)+anisok(n,i,j,k)*wei  
            endfor
            if (iopt eq 1 and iunit eq 1) $    ; differential flux
                then plsfluxa(n,i,j)=plsfluxa(n,i,j)/weiSum
            if (iopt eq 2) then plsfluxa(n,i,j)=pf*plsfluxa(n,i,j)
            if (iopt eq 3) then plsfluxa(n,i,j)=nf*plsfluxa(n,i,j)
            if (iopt eq 4 and denSum gt 0.) then plsfluxa(n,i,j)= $
                                            Tf*plsfluxa(n,i,j)/denSum ; T in keV
            anisoa(n,i,j)=anisoa(n,i,j)/weiSum
            if (ilog eq 1) then begin
               if (plsfluxa(n,i,j) eq 0.) then plsfluxa(n,i,j)=-50.
               if (plsfluxa(n,i,j) gt 0.) then $
                   plsfluxa(n,i,j)=alog10(plsfluxa(n,i,j))
            endif
        endfor
    endfor
endfor

; setup plot ranges
yon=' '
fmax=max(plsfluxa)    
print,' fmax = ',fmax
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmax => ',fmax
fmin=0.   
if (ilog eq 1) then fmin=fmax-3.   
print,' fmin = ',fmin
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmin => ',fmin
amax=1.
amin=-1.

; smooth data
read,' Do you want to smooth the data? (y/n) => ',yon
nframe=0
if (yon eq 'y') then read,' enter number of intermediate frame => ',nframe
ntnew=(ntime-1)*(nframe+1)+1
ronew=fltarr(ntnew,nr,nmlt+1)  &  mltonew=ronew  
anisonew=ronew                 &  plsfluxnew=ronew
if (nframe eq 0) then begin
   hournew=houra        
   ronew=roa        
   mltonew=mltoa   
   anisonew=anisoa     
   plsfluxnew=plsfluxa
endif else begin
   hournew=interpol(houra,ntnew)
   for i=0,nr-1 do begin 
       for j=0,nmlt do begin
           ronew(*,i,j)=interpol(roa(*,i,j),ntnew)
           mltonew(*,i,j)=interpol(mltoa(*,i,j),ntnew)
           anisonew(*,i,j)=interpol(anisoa(*,i,j),ntnew)
           plsfluxnew(*,i,j)=interpol(plsfluxa(*,i,j),ntnew)
       endfor
   endfor
endelse

; Setup window
ips=0
set_plot,'x'
x_wsize=684
y_wsize=480
window,1,xpos=200,ypos=100,xsize=x_wsize,ysize=y_wsize

; Plot plasma flux (or pressure, or density) and pitch angle anisotropy
for n=0,ntnew-1 do begin
    plot_flux:
    polyfill,[0.,1.,1.,0.,0.],[0.,0.,1.,1.,0.],color=255,/normal  ; white BG
    hour=hournew(n)
    ro(*,*)=ronew(n,*,*)
    mlto(*,*)=mltonew(n,*,*)
    aniso(*,*)=anisonew(n,*,*)
    plsflux(*,*)=plsfluxnew(n,*,*)

    ; plot anisotropy
    xi=0.56
    if (min(plsflux) lt fmin) then aniso(where(plsflux lt fmin))=0.
    bar_lab='pitch angle anisotropy'
    if (iopt eq 2) then bar_lab='pressure anisotropy'
    if (iopt eq 4) then bar_lab='temperature anisotropy'
    plot_equator,aniso,amax,amin,xi,bar_lab,halfl,ro,mlto,lat,x_wsize,y_wsize, $
                 dir,blackc,g_thick,iround

    ; plot plasma fluxes or ion pressure, or density
    xi=0.07
    if (iopt eq 1) then begin 
       if (iunit eq 1) then bar_lab='flux (/keV/cm2/sr/s)'
       if (iunit eq 2) then bar_lab='flux (/cm2/sr/s)'
    endif
    if (iopt eq 2) then bar_lab='pressure (nPa)'
    if (iopt eq 3) then bar_lab='density (cm^-3)'
    if (iopt eq 4) then bar_lab='temperature (keV)'
    if (ilog eq 1) then bar_lab='log '+bar_lab
    plot_equator,plsflux,fmax,fmin,xi,bar_lab,halfl,ro,mlto,lat,x_wsize, $
                 y_wsize,dir,blackc,g_thick,iround

    ;label the plot
    hours=string(fix(hour),'(i3)')
    minute=round((hour-float(fix(hour)))*60.)
    minutes=string(minute,'(i2.2)')
    elabel=string(e0,flab)+' - '+string(e1,flab)
    xyouts,0.07,0.93,fhead,size=2.,color=blackc,/normal
    xyouts,0.5,0.93,hours+':'+minutes+' UT',size=2.,alignment=0.5, $
           color=blackc,/normal 
    xyouts,0.5,0.85,elabel+' keV '+species,size=2.,alignment=0.5, $
           color=blackc,/normal

    ; Make gif file
    if (ips eq 1) then goto,close_gif
    gimg = tvrd(0,0)
    write_gif,fhead+'_'+fmiddle+'.gif',gimg,red(0:255),green(0:255), $ 
              blue(0:255),/multiple

    ; Make ps file if (n eq ntnew-1)
    if (n eq ntnew-1) then begin
       set_plot,'ps'
       ips=1
       !p.charthick=3.0
       g_thick=3.0
       device,filename='plot_fls.ps',/inches,/color,yoffset=0.5,xoffset=0.3, $
              xsize=5.*x_wsize/y_wsize,ysize=5.,bits_per_pixel=24 
       goto,plot_flux
    endif
endfor

close_gif:
write_gif,fhead+'_'+fmiddle+'.gif',/close
device,/close
read,'Do you want to continue?  (y)yes, (n)no => ',yon
if (yon eq 'y') then goto,new_plot

end
