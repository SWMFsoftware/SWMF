;-------------------------------------------------------------------------------
;
;                              plot_fls_L_Time.pro
;
; IDL procedure reads *.fls file and plot the pitch-angle local-time averaged
; flux in a L-time plot.
;
; Created on 9 November 2005 by Mei-Ching Fok, Code 612.2, NASA GSFC.
;-------------------------------------------------------------------------------

close,/all
;write_gif,/close

; get colors for color table palette
loadct,39       ;  rainbow with white
tvlct, red, green, blue, /get

; read file name and energy, sina, lat information 
fhead=' '
read,'enter the file head (e.g., 2000_225_h) => ',fhead
openr,2,fhead+'.fls'
readf,2,rc,ir,ip,je,ig   ; je=number of energy grid, ig=number of pa grid
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
energy=fltarr(je)    &   sina=fltarr(ig)   &   lat=fltarr(ir)   &   sinao=sina
Ebound=fltarr(je+1)
readf,2,energy
;readf,2,Ebound
readf,2,sina
readf,2,lat
read,'plot flux at  1 = equator,  2 = ionosphere => ',iplt
y_title='L'
if (iplt eq 1) then begin
   read,'plot at?  1=L-invariant, 2=L-geocentric, 4=L* => ',iLs
   if (iLs eq 1) then y_title='L-invariant' 
   if (iLs eq 2) then y_title='L geocentric'
   if (iLs eq 4) then y_title='L*'  
endif
alt=9999.
if (iplt eq 2) then read,'altitude (km) => ',alt
ri=(alt+6375.)/6375.           ; ionosphere distance in RE
rim=ri*6.375e6                 ; ionosphere distance in m
xme=7.9e15                     ; magnetic dipole moment of the Earth
read,'number of data set in time => ',ntime
read,'flux in?  0 = linear scale, or 1 = log scale => ',ilog

; Calculate L-values
Lshell=lat
for i=0,ir-1 do begin
    cosl=cos(lat(i)*!pi/180.)
    Lshell(i)=rc/cosl/cosl
endfor

; define energy labels
for k=1,je-1 do Ebound(k)=sqrt(energy(k-1)*energy(k))
Ebound(0)=energy(0)*energy(0)/Ebound(1)
Ebound(je)=energy(je-1)*energy(je-1)/Ebound(je-1)
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
plsflux3=fltarr(nr,nmlt,je) ;& plsflux2=fltarr(nmlt,je)
plsfluxk=fltarr(ntime,nr,je)    
plsflux=fltarr(ntime,nr)      
hour=fltarr(ntime)
fy=fltarr(ig)      &    fynew=fy    &     dmu=fy        &    cosa=fy
cosa(*)=cos(asin(sina(*)))
if (iLs eq 2) then ro=fltarr(nr,nmlt)
for m=0,ig-1 do begin
    if (m eq 0) then sina0=0.
    if (m gt 0) then sina0=0.5*(sina(m)+sina(m-1))
    if (m eq (ig-1)) then sina1=1.
    if (m lt (ig-1)) then sina1=0.5*(sina(m)+sina(m+1))
    dmu(m)=sqrt(1.-sina0*sina0)-sqrt(1.-sina1*sina1)
    ;print,  dmu(m),sina0,sina0,sina1,sina1
endfor
; setup arrays for L*
if (iLs eq 4) then begin
   Lstar=fltarr(nr,nmlt) & mlto=Lstar
   Lstar_max=fltarr(ntime) & Lstar1=fltarr(nr,nmlt+1)
   xe=fltarr(nr,nmlt+1) & ye=xe
; Half-shaded Earthi
   x_1=fltarr(61) & y_1=x_1 & earth_x=fltarr(62) & earth_y=earth_x
   for i=-30,30 do begin
       earth_x(i+30)=cos(3*i*!dtor)
       earth_y(i+30)=sin(3*i*!dtor)
       x_1(i+30)=cos(i*6.*!dtor)
       y_1(i+30)=sin(i*6.*!dtor)
   endfor
; setup levels
   llevel=9
   llvl=findgen(llevel)  &  lcolr=intarr(llevel)
   colrmax=254           
   colrmin=1
   ncolor=colrmax-colrmin+1
   dcolr=(float(colrmax)-float(colrmin))/(llevel-1)
   for i=0,llevel-1 do lcolr(i)=round(float(colrmin)+i*dcolr)
   lcolr(4)=lcolr(3)
   lcolr(3)=round(0.5*(float(lcolr(2)+lcolr(3))))
; position to plot L*
   x_wsize=570
   y_wsize=446
   pos1=[0.12,0.08,0.76,0.92]
   tickn=['10','5','0','-5','-10']
   tickn0=[' ',' ',' ',' ',' ']
   set_plot,'ps'
   file_name=fhead+'_Lstar.ps'
   device,filename=file_name,/inches,/color,yoffset=3.0, $
          xoffset=0.3,xsize=7.,ysize=7.*y_wsize/x_wsize
endif ; end if iLs eq 4 loop
dE=fltarr(je)
for k=0,je-1 do dE(k)=Ebound(k+1)-Ebound(k)

; Read fluxes and calculate pa-mlt averaged flux
Lmax0=15.
for n=0,ntime-1 do begin
    readf,2,hour1,Lmax0
    if (iLs eq 4) then Lstar_max(n)=Lmax0
    hour(n)=hour1
    print,n,hour1,Lmax0
    for i=0,nr-1 do begin
        bi=sqrt(4.-3.*ri/Lshell(i))*xme/rim^3  ; assume dipole
        for j=0,nmlt-1 do begin
            readf,2,lat1,xmlt,ro0,xmlto,bo,Lstar0    ; bo: equatorial B in T
            if (iLs eq 2) then ro(i,j)=ro0
            if (iLs eq 4) then begin
	       if (Lstar0 gt 0.) then Lstar(i,j)=Lstar0 else Lstar(i,j)=-1.*Lstar0
               mlto(i,j)=xmlto
               xe(i,j)=ro0*cos(xmlto*15.*!dtor)
               ye(i,j)=ro0*sin(xmlto*15.*!dtor)
            endif
            sqrtbobi=sqrt(bo/bi)
            sinao(*)=sina(*)*sqrtbobi
            for k=0,je-1 do begin
                readf,2,fy
                fynew=fy
                if (iplt eq 2) then fynew=interpol(fy,sina,sinao)
                plsflux3(i,j,k)=0.    ; pitch-angle averaged flux
                for m=0,ig-1 do begin
                    plsflux3(i,j,k)=plsflux3(i,j,k)+fynew(m)*dmu(m)
                    ;print,plsflux2(j,k),fynew(m),dmu(m)
                endfor
            endfor
        endfor                          ; plasfluxk: pa-mlt averaged flux
        ;for k=0,je-1 do begin
        ;    plsfluxk(n,i,k)=total(plsflux2(*,k))/nmlt 
        ;endfor
        if (iLs eq 4) then begin
           xe(i,nmlt)=xe(i,0) & ye(i,nmlt)=ye(i,0)
        endif
    endfor
    ; Draw L* contours
    if (iLs eq 4) then begin            
       Lstar1(*,0:nmlt-1)=Lstar(*,0:nmlt-1) & Lstar1(*,nmlt)=Lstar1(*,0)
       dmax=fix(Lstar_max(n))
       if fix(Lstar_max(n)) lt llevel-1 then begin
          llvl(fix(Lstar_max(n))+1)=Lstar_max(n)
          lcolr(fix(Lstar_max(n))+1:llevel-1)=255
       endif
       clvl=findgen(dmax-1)+2.
       Lmax_lab=string(Lstar_max(n),'(f5.2)')
       min1=fix((hour1-fix(hour1))*60)
       UT_lab=string(hour1,'(i3.3)')+':'+string(min1,'(i2.2)')
       xyouts,0.12,0.96,UT_lab+' hh:mm',size=2.0,color=0,/normal
       contour,Lstar1,xe,ye,xticks=4,yticks=4,xtickn=tickn,ytickn=tickn, $
        levels=Llvl,pos=pos1,charsize=2,xrange=[-10,10],yrange=[-10,10],xstyle=1,ystyle=1, $
        c_colors=lcolr,xtitle='RE',background=255,color=0,/fill,/noerase
       contour,Lstar1,xe,ye,xticks=4,yticks=4,xtickn=tickn,ytickn=tickn, $
        levels=clvl,c_labels=0,pos=pos1,charsize=2,c_thick=2, $
        xrange=[-10,10],yrange=[-10,10],xstyle=1,ystyle=1, $
        c_colors=0,xtitle='RE',background=255,color=0,/noerase
       contour,Lstar1,xe,ye,xticks=4,yticks=4,xtickn=tickn,ytickn=tickn, $
        levels=Lstar_max(n),c_labels=0,pos=pos1,charsize=2,c_thick=4, $
        xrange=[-10,10],yrange=[-10,10],xstyle=1,ystyle=1, $
        c_colors=0,xtitle='RE',background=255,color=0,/noerase
       axis,0,0,xaxis=1,xrange=[10,-10],xticks=ntick,color=0,xstyle=1,xtickn=tickn0
       axis,0,0,yaxis=1,yrange=[10,-10],yticks=ntick,color=0,ystyle=1,ytickn=tickn0
       oplot,xe(nr-1,0:nmlt-1),ye(nr-1,0:nmlt-1)
       oplot,x_1,y_1
       polyfill,earth_x,earth_y
       xyouts,2.,Lstar_max(n)+0.5,'L* max = '+Lmax_lab,charsize=1.75
     ; draw color bar and label
       for i=fix(Lstar_max(n)+1),llevel-1 do lcolr(i)=round(float(colrmin)+i*dcolr)
       llvl=findgen(llevel)
       x1=0.76+0.14
       x2=x1+0.03
       dy=0.8/(llevel-1)
       y2=0.08
       for i=1,llevel-1 do begin
           y1=y2
           y2=y1+dy*1.01
           polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=lcolr(i),/normal  
           xyouts,x2,y1-0.02,string(llvl(i),'(i1.1)'),size=2.0,color=0,/normal
       endfor
       xyouts,x2,0.88,string(llvl(llevel-1)+1,'(i1.1)'),size=2.0,color=0,/normal
       xyouts,x2+0.04,0.5,'L*',alignment=0.5,orientation=90.,color=0,size=1.5, $
           /normal
       erase
    endif   ; end of if iLs eq 4 loop
    ; map flux to fixed L grid, geometric L, or L* 
    for j=0,nmlt-1 do begin
        if (iLs eq 2) then ii=interpol(findgen(nr),ro(*,j),Lshell)
        if (iLs eq 4) then ii=interpol(findgen(nr),Lstar(*,j),Lshell)
        if (iLs gt 1) then begin
           for k=0,je-1 do begin
               plsflux1=interpolate(plsflux3(*,j,k),ii,missing=0.)
               plsflux3(*,j,k)=plsflux1(*)
           endfor
        endif 
    endfor
    ; find pa-mlt averaged flux
    for i=0,nr-1 do begin
        plsfluxk(n,i,0:je-1)=0.   
        if (Lshell(i) gt ri or iplt eq 1) then begin
           for k=0,je-1 do plsfluxk(n,i,k)=total(plsflux3(i,*,k))/nmlt 
           if (iLs eq 4 and Lshell(i) gt Lmax0) then plsfluxk(n,i,*)=0.
        endif
    endfor
endfor
close,2
if (iLs eq 4) then device,/close_file

new_plot:
; choose which energy to be displaced
print,'energy bin (keV): '
for i=0,je-1 do print,i,Ebound(i:i+1),format='("  (",i2,") ",f6.1," - ",f6.1)'
read, 'lower and upper energy bins (ie1,ie2) => ',ie1,ie2
read,'unit? (1)flux/(keV cm2 s sr), (2)flux/(cm2 s sr) => ',iunit
e0=Ebound(ie1)
e1=Ebound(ie2+1)
fmiddle=string(round(e0),ilab)+'-'+string(round(e1),ilab)+'keV'

; calculate total flux 
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        plsflux(n,i)=0.
        weiSum=0.
        for k=ie1,ie2 do begin    
            weiSum=weiSum+dE(k)
            ;print,plsfluxk(n,i,k),dE(k)
            plsflux(n,i)=plsflux(n,i)+plsfluxk(n,i,k)*dE(k)
        endfor
        ;print,plsflux(n,i),weiSum
        if (iunit eq 1) then plsflux(n,i)=plsflux(n,i)/weiSum  ; diff flux
    endfor
endfor
if (ilog eq 1) then plsflux=alog10(plsflux)

; setup plot ranges
yon=' '
fmax=max(plsflux)    
print,' fmax = ',fmax
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmax => ',fmax
fmin=0.   
if (ilog eq 1) then fmin=fmax-5.
if (min(plsflux) lt fmin) then plsflux(where(plsflux lt fmin))=fmin

; setup levels
nlevel=59
lvl=fltarr(nlevel)    &    colr=intarr(nlevel)
dlvl=(fmax-fmin)/(nlevel-1)
colrmax=254           
colrmin=1
ncolor=colrmax-colrmin+1
dcolr=(float(colrmax)-float(colrmin))/(nlevel-1)
for i=0,nlevel-1 do begin
    lvl(i)=fmin+i*dlvl
    colr(i)=round(float(colrmin)+i*dcolr)
endfor

; Setup window
ips=0
set_plot,'x'
x_wsize=570
y_wsize=446
window,1,xpos=200,ypos=100,xsize=x_wsize,ysize=y_wsize
!p.color=0
!p.background=255

; Plot plasma flux 
x0=0.1
xf=0.8
y0=0.13 
yf=0.83
ym=0.5*(y0+yf)
elabel=string(e0,flab)+' - '+string(e1,flab)
mtitle=elabel+' keV '+species
read,'Do you want to reduce the temporal resolution? y(yes), n(no) => ',yon
if (yon eq 'n') then begin
   plsfluxs=plsflux
   hours=hour
endif else begin
   read,'enter number of time sector => ',ntsec
   ntnew=ntsec*2
   dn=fix(ntime/ntsec)
   plsfluxs=fltarr(ntnew,nr)   &   hours=fltarr(ntnew)
   for n=0,ntsec-1 do begin
       n0=n*dn
       n1=n0+dn
       n2=n*2
       hours(n2)=hour(n0)   
       hours(n2+1)=hour(n1)-0.01
       for i=0,nr-1 do plsfluxs(n2,i)=(total(plsflux(n0:n1,i)))/(n1-n0+1)
       plsfluxs(n2+1,0:nr-1)=plsfluxs(n2,0:nr-1)  
   endfor
endelse
plotFlux:
x_tick=fix(max(hour))/24
contour,plsfluxs,hours,Lshell,yrange=[1,8],ystyle=1,xtitle='hour',ytitle=y_title, $
        xrange=[min(hour),max(hour)],xstyle=1,xticks=x_tick,title=mtitle,  $
        charsize=1.5,pos=[x0,y0,xf,yf],levels=lvl,c_colors=colr,/fill 
xyouts,x0,0.93,strmid(fhead,0,8),size=2.,/normal

; draw color bar and label
x1=xf+0.06
x2=x1+0.03
dy=(yf-y0)/(colrmax-colrmin+1)
for i=colrmin,colrmax do begin
    y1=y0+i*dy
    y2=y1+dy*1.02
    polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=i,/normal  
endfor
if (iunit eq 1) then bar_lab='flux (/keV/cm2/sr/s)'
if (iunit eq 2) then bar_lab='flux (/cm2/sr/s)'
if (ilog eq 1) then bar_lab='log '+bar_lab
xyouts,x2,y0,string(fmin,'(f4.1)'),size=1.5,/normal
xyouts,x2,yf,string(fmax,'(f4.1)'),size=1.5,/normal
xyouts,x2+0.03,ym,bar_lab,alignment=0.5,orientation=90.,size=1.5,/normal
if (iplt eq 1) then iplt_lab='equatorial flux'
if (iplt eq 2) then iplt_lab=string(fix(alt),format='(i4)')+' km altitude'
xyouts,0.5*(x0+xf),yf+0.08,iplt_lab,alignment=0.5,size=1.1,/normal

; Make gif file
if (ips eq 1) then goto,newPlot
gimg = tvrd(0,0)
;write_gif,fhead+'_'+fmiddle+'.gif',gimg,red(0:255),green(0:255),blue(0:255) 

; Make ps file
set_plot,'ps'
ips=1
device,filename='plot_fls_L_Time.ps',/inches,/color,yoffset=0.5,xoffset=0.3, $
       xsize=7.,ysize=7.*y_wsize/x_wsize
goto,plotFlux

newPlot:
read,'Do you want to continue?  (y)yes, (n)no => ',yon
if (yon eq 'y') then goto,new_plot
;write_gif,fhead+'_'+fmiddle+'.gif',/close
device,/close

end
