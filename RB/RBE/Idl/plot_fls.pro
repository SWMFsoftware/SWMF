;-------------------------------------------------------------------------------
;
;                                plot_fls.pro
;
; IDL procedure reads *.fls file and plot the pitch-angle averaged flux and 
; the pitch-angle distribution at the equatorial plane.  This program can
; handle multiple data sets in time.
;
; Created on 10 September 2002 by Mei-Ching Fok, NASA GSFC LEP, Code 692
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
pro plot_equator,var,vmax,vmin,xi,bar_lab,halfl,ro,mlto,x_wsize,y_wsize,dir
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

; force var equal or gt vmin for flux and pressure
if (min(var) lt vmin) then var(where(var le vmin))=vmin

; polyfill background color in black
polyfill,[xi,xf,xf,xi,xi],[yi,yi,yf,yf,yi],/normal

; plot flux or PA anisotropy
if (dir eq 1) then phi=mlto               ; Sun to the left
if (dir eq 2) then phi=mlto+!pi           ; Sun to the right
polar_contour,var,phi,ro,xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
      levels=lvl,c_colors=colr,xstyle=1,ystyle=1,xticks=4,yticks=4, $
      xtitle='RE',pos=[xi,yi,xf,yf],/fill,/noerase

; Add color bar and label
xb=xi
yb=0.090
x0b=xb*x_wsize
y0b=yb*y_wsize
temp=bytarr(ncolor,2)
for i=0,ncolor-1 do begin
   temp(i,0)=i+colrmin
   temp(i,1)=i+colrmin
endfor
bar_len=200
color_bar=congrid(temp,bar_len,20)
tv,color_bar,x0b,y0b
v_min=string(vmin,'(f4.1)')
if (vmax lt 100.) then v_max=string(vmax,'(f4.1)')
if (vmax ge 100.) then v_max=string(vmax,'(e7.1)')
if (bar_lab eq 'pitch angle anisotropy') then v_min='field aligned'
if (bar_lab eq 'pitch angle anisotropy') then v_max='perpendicular'
xyouts,xb,0.6*yb,v_min,alignment=0.5,size=1.2,/normal
xlab=xb+float(bar_len)/x_wsize
xyouts,xlab,0.6*yb,v_max,alignment=0.5,size=1.2,/normal
xyouts,0.5*(xlab+xb),1.7*yb,bar_lab,alignment=0.5,size=1.1,/normal

; fill cirle of radius ro(0,0) centered the earth with color2
color2=colrmin
if (bar_lab eq 'pitch angle anisotropy') then begin
   m=(nlevel-1)/2
   color2=colr(m)
endif

; Draw earth , geosynchronous
npt=400                             ; no. of points in a circle
npt2=npt/2
night_x=fltarr(npt2+2)     &    day_x=night_x
night_y=fltarr(npt2+2)     &    day_y=night_y
e1x = fltarr(npt+1)
e1y = fltarr(npt+1)
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
oplot,e1x,e1y,color=255,thick=1.5         ; plot geosynchronous
polyfill,ro(0,0)*night_x,ro(0,0)*night_y,color=color2
polyfill,ro(0,0)*day_x,ro(0,0)*day_y,color=color2
polyfill,night_x,night_y
polyfill,day_x,day_y,color=255

return
end


;-----------------------------------------------------------------------------
; main routine
;-----------------------------------------------------------------------------

close,/all

; get colors for color table palette
read,'colar table, 1=rainbow with white, 2=HENA-type => ',icr
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

; read file name and energy, sina, lat information 
fhead=' '
read,'enter the file head (e.g., 2000_225_h) => ',fhead
openr,2,fhead+'.fls'
readf,2,rc,ir,ip,je,ig   ; je=number of energy grid, ig=number of pa grid
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
energy=fltarr(je)    &   sina=fltarr(ig)   &   lat=fltarr(ir)
Ebound=fltarr(je+1)
readf,2,energy
readf,2,Ebound
readf,2,sina
readf,2,lat
read,'number of data set in time => ',ntime
read,'flux or pressure in?  0 = linear scale, or 1 = log scale => ',ilog
read,'Sun to the?   1 = left,  2 = right => ',dir
read,'half length of the plot => ', halfl

; define energy labels
if (Ebound(je) lt 999.5) then ilab='(i3.3)'
if (Ebound(je) ge 999.5) then ilab='(i4.4)'
if (Ebound(je) lt 999.5) then flab='(f5.1)'
if (Ebound(je) ge 999.5) then flab='(f6.1)'
 
; setup energy and species
spe=strmid(fhead,9,1)                ; species
if (spe eq 'h') then species='H+'
if (spe eq 'o') then species='O+'
if (spe eq 'e') then species='e-'

; setup arrays
plsfluxk=fltarr(ntime,nr,nmlt+1,je)    &  anisok=plsfluxk
plsfluxa=fltarr(ntime,nr,nmlt+1)       &  anisoa=plsfluxa
roa=fltarr(ntime,nr,nmlt+1)   &  mltoa=roa   
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

; factor of calculating ion pressure
if (spe eq 'h') then mass=1.
if (spe eq 'o') then mass=16.
if (spe eq 'e') then mass=5.4462e-4
pf=4.*!pi*1.e4/3.*sqrt(2.*1.673e-27*mass)*sqrt(1.6e-16)*1.e9  ; pressure in nPa

; Read fluxes and calculate PA anisotropy
for n=0,ntime-1 do begin
    readf,2,hour,format='(8x,f6.2)'
    houra(n)=hour
    print,n,hour
    for i=0,nr-1 do begin
        for j=0,nmlt-1 do begin
            readf,2,lat1,mlt,ro1,mlto1
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
; choose which energy to be displaced
print,'energy bin (keV): '
for i=0,je-1 do print,i,Ebound(i:i+1),format='("  (",i2,") ",f6.1," - ",f6.1)'
read, 'lower and upper energy bins (ie1,ie2) => ',ie1,ie2
read,'plot? (1)flux,  (2)pressure => ',iopt
if (iopt eq 1) then $
    read,'unit? (1)flux/(keV cm2 s sr), (2)flux/(cm2 s sr) => ',iunit
e0=Ebound(ie1)
e1=Ebound(ie2+1)
fmiddle=string(round(e0),ilab)+'-'+string(round(e1),ilab)+'keV'
if (iopt eq 1) then fmiddle=fmiddle+'_flux' 
if (iopt eq 2) then fmiddle=fmiddle+'_pressure'

; calculate total flux or pressure
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        for j=0,nmlt do begin

            plsfluxa(n,i,j)=0.
            anisoa(n,i,j)=0.
            weiSum=0.
            for k=ie1,ie2 do begin    
                if (iopt eq 1) then wei=dE(k)
                if (iopt eq 2) then wei=dE(k)*sqrt(energy(k))
                weiSum=weiSum+wei
                plsfluxa(n,i,j)=plsfluxa(n,i,j)+plsfluxk(n,i,j,k)*wei  
                anisoa(n,i,j)=anisoa(n,i,j)+anisok(n,i,j,k)*wei  
            endfor
            if (iopt eq 1 and iunit eq 1) $    ; differential flux
                then plsfluxa(n,i,j)=plsfluxa(n,i,j)/weiSum
            if (iopt eq 2) then plsfluxa(n,i,j)=pf*plsfluxa(n,i,j)
            anisoa(n,i,j)=anisoa(n,i,j)/weiSum
        endfor
    endfor
endfor
if (ilog eq 1) then plsfluxa=alog10(plsfluxa)

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

; Setup window
!p.charsize=3
x_wsize=1140
y_wsize=800
window,1,xpos=200,ypos=100,xsize=x_wsize,ysize=y_wsize
!p.color=blackc

; Plot plasma flux (or pressure) and pitch angle anisotropy
for n=0,ntime-1 do begin
    polyfill,[0.,1.,1.,0.,0.],[0.,0.,1.,1.,0.],color=255,/normal  ; white BG
    hour=houra(n)
    ro(*,*)=roa(n,*,*)
    mlto(*,*)=mltoa(n,*,*)
    aniso(*,*)=anisoa(n,*,*)
    plsflux(*,*)=plsfluxa(n,*,*)

    ; plot PA anisotropy
    xi=0.56
    if (min(plsflux) lt fmin) then aniso(where(plsflux lt fmin))=0.
    bar_lab='pitch angle anisotropy'
    plot_equator,aniso,amax,amin,xi,bar_lab,halfl,ro,mlto,x_wsize,y_wsize,dir

    ; plot plasma fluxes or ion pressure
    xi=0.07
    if (iopt eq 1) then begin 
       if (iunit eq 1) then bar_lab='flux (/keV/cm2/sr/s)'
       if (iunit eq 2) then bar_lab='flux (/cm2/sr/s)'
    endif
    if (iopt eq 2) then bar_lab='pressure (nPa)'
    if (ilog eq 1) then bar_lab='log '+bar_lab
    plot_equator,plsflux,fmax,fmin,xi,bar_lab,halfl,ro,mlto,x_wsize,y_wsize,dir

    ;label the plot
    hours=string(fix(hour),'(i3)')
    minute=round((hour-float(fix(hour)))*60.)
    minutes=string(minute,'(i2.2)')
    elabel=string(e0,flab)+' - '+string(e1,flab)
    xyouts,0.07,0.93,fhead,size=2.,/normal
    xyouts,0.5,0.93,hours+':'+minutes+' UT',size=2.,alignment=0.5,/normal 
    xyouts,0.5,0.85,elabel+' keV '+species,size=2.,alignment=0.5,/normal

    ; Make gif file
    gimg = tvrd(0,0,order=0,true=1)
    
    ;write_tiff,fhead+'_'+fmiddle+strtrim(n,2)+'.tif',gimg,red=red(0:255),$
    ;  green=green(0:255),blue=blue(0:255) ;,/multiple

    write_image,fhead+'_'+fmiddle+strtrim(n,2)+'.tif','TIFF',gimg,red=red(0:255),$
      green=green(0:255),blue=blue(0:255) ;,/multiple
endfor

;write_tiff,fhead+'_'+fmiddle+'.tiff',red=red(0:255),$
;      green=green(0:255),blue=blue(0:255);,/close
read,'Do you want to continue?  (y)yes, (n)no => ',yon
if (yon eq 'y') then goto,new_plot

end
