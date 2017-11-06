;-------------------------------------------------------------------------------
;
;                              plot_fls_L_Time.pro
;
; IDL procedure reads *.fls file and plot the pitch-angle local-time averaged
; psd in a L-time plot.
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
openr,2,fhead+'.psd'
readf,2,rc,ir,ip,nm,nk   ; je=number of energy grid, ig=number of pa grid
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
xk=fltarr(nk) & xmm=fltarr(nm) & lat=fltarr(ir)
readf,2,xk
;; readf,2,Ebound
readf,2,xmm
readf,2,lat
;; read,'plot psd at  1 = equator,  2 = ionosphere => ',iplt
iplt=1
alt=9999.
if (iplt eq 2) then read,'altitude (km) => ',alt
ri=(alt+6375.)/6375.           ; ionosphere distance in RE
rim=ri*6.375e6                 ; ionosphere distance in m
xme=7.9e15                     ; magnetic dipole moment of the Earth
read,'number of data set in time => ',ntime
read,'psd in?  0 = linear scale, or 1 = log scale => ',ilog

; Calculate L-values
Lshell=lat
for i=0,ir-1 do begin
    cosl=cos(lat(i)*!pi/180.)
    Lshell(i)=rc/cosl/cosl
endfor

; define energy labels
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
psd=fltarr(ntime,nr,nmlt,nm,nk)
plspsd2=fltarr(nmlt,nm,nk)
plspsdk=fltarr(ntime,nr,nm,nk)
plspsd=fltarr(ntime,nr)
roa=fltarr(ntime,nr,nmlt+1)   &  mltoa=roa   
;irm=intarr(ntime,nmlt)        &  iba=irm
ro=fltarr(nr,nmlt+1)          &  mlto=ro
hour=fltarr(ntime)
psdy=fltarr(nk)
;; fy=fltarr(ig)      &    fynew=fy    &     dmu=fy        &    cosa=fy
;; cosa(*)=cos(asin(sina(*)))
;; for m=0,ig-1 do begin
;;     if (m eq 0) then sina0=0.
;;     if (m gt 0) then sina0=0.5*(sina(m)+sina(m-1))
;;     if (m eq (ig-1)) then sina1=1.
;;     if (m lt (ig-1)) then sina1=0.5*(sina(m)+sina(m+1))
;;     dmu(m)=sqrt(1.-sina0*sina0)-sqrt(1.-sina1*sina1)
;;     ;print,  dmu(m),sina0,sina0,sina1,sina1
;; endfor
;; dE=fltarr(je)
;; for k=0,je-1 do dE(k)=Ebound(k+1)-Ebound(k)

; Read psd and calculate mlt averaged psd
for n=0,ntime-1 do begin
   readf,2,hour1
   hour(n)=hour1
   print,n,hour1
   for i=0,nr-1 do begin
      ;; bi=sqrt(4.-3.*ri/Lshell(i))*xme/rim^3 ; assume dipole
      for j=0,nmlt-1 do begin
         readf,2,lat1,mlt,ro1,mlto1,bo ;,irm1,iba1
         ;; sqrtbobi=sqrt(bo/bi)
         ;; sinao(*)=sina(*)*sqrtbobi
         for m=0,nm-1 do begin
            readf,2,psdy
            psd(n,i,j,m,0:nk-1)=psdy
            ;; psdynew=psdy
            ;; if (iplt eq 2) then fynew=interpol(fy,sina,sinao)
            ;; plspsd2(j,k)=0.    ; pitch-angle averaged psd
            ;; for k=0,nk-1 do begin
            ;;    plspsdk(j,m,0:nk-1)=psdy
            ;; plspsd2(j,m,k)=plspsd2(j,k)+fynew(m)*dmu(m)
            ;; print,plspsd2(j,k),fynew(m),dmu(m)
            ;; endfor
         endfor
      endfor

; plspsdk: mlt averaged psd
      for m=0,nm-1 do begin
         for k=0,nk-1 do begin
            plspsdk(n,i,m,k)=total(psd(n,i,*,m,k))/nmlt
         endfor
      endfor
   endfor
endfor

close,2

new_plot:
; choose which mu and K to be displayed
print,'mu bins (nT^0.5 R_E): '
for m=0,nm-1 do print,m,xmm(m),format='("  (",i2,") ",e8.2," keV/nT")'
read, 'select mu bin (im1) => ',im1
mu0=xmm(im1)
if (mu0 lt 0000.01) then ilab='(e7.1)'
if (mu0 ge 0000.01) then ilab='(f4.2)'
if (mu0 ge 0000.10) then ilab='(f3.2)'
if (mu0 ge 0001.00) then ilab='(f3.1)'
if (mu0 ge 0010.00) then ilab='(f4.1)'
if (mu0 ge 0100.00) then ilab='(f5.1)'
if (mu0 ge 1000.00) then ilab='(f6.1)'
if (mu0 lt 0000.01) then flab='(e7.1)'
if (mu0 ge 0000.01) then flab='(f4.2)'
if (mu0 ge 0000.10) then flab='(f3.2)'
if (mu0 ge 0001.00) then flab='(f3.1)'
if (mu0 ge 0010.00) then flab='(f4.1)'
if (mu0 ge 0100.00) then flab='(f5.1)'
if (mu0 ge 1000.00) then flab='(f6.1)'
print,'K bins (nT^0.5 R_E): '
for k=0,nk-1 do print,k,xk(k),format='("  (",i2,") ",f7.1," nT^0.5 R_E")'
read, 'select K bin (ik1) => ',ik1
k0=xk(ik1)
if (k0 lt 00010.) then fklab='(f3.1)'
if (k0 ge 00010.) then fklab='(f4.1)'
if (k0 ge 00100.) then fklab='(f5.1)'
if (k0 ge 01000.) then fklab='(f6.1)'
if (k0 ge 10000.) then fklab='(f7.1)'
print,"mu0: ",mu0
print,"K0: ",K0
plspsd=plspsdk(*,*,im1,ik1)
iunit=0
read,'psd units? (1) c3/(MeV3 cm3), (2) s3/(kg3 m6) => ',iunit
if (iunit eq 2) then plspsd=plspsd/1.e18/q_e^3*c_speed^3
fmiddle=string(round(mu0),ilab)+'keV/nT'+string(k0,fklab)+'nT^{0.5} R_E'

; calculate total psd 
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        ;; plspsd(n,i)=0.
        ;; weiSum=0.
        ;; for k=ie1,ie2 do begin    
        ;;     weiSum=weiSum+dE(k)
        ;;     ;print,plspsdk(n,i,k),dE(k)
        ;;     plspsd(n,i)=plspsd(n,i)+plspsdk(n,i,k)*dE(k)
        ;; endfor
        ;; ;print,plspsd(n,i),weiSum
        ;; if (iunit eq 1) then plspsd(n,i)=plspsd(n,i)/weiSum  ; diff psd
       if (ilog eq 1) then begin
          if (plspsd(n,i) eq 0.) then plspsd(n,i)=-50.
          if (plspsd(n,i) gt 0.) then $
             plspsd(n,i)=alog10(plspsd(n,i))
       endif
    endfor
endfor
;if (ilog eq 1) then plspsd=alog10(plspsd)

; setup plot ranges
yon=' '
fmax=max(plspsd)    
print,' fmax = ',fmax
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmax => ',fmax
fmin=0.   
if (ilog eq 1) then fmin=fmax-6.
print,' fmin = ',fmin
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmin => ',fmin
;if (min(plspsd) lt fmin) then plspsd(where(plspsd lt fmin))=fmin

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

; Plot plasma psd 
x0=0.1
xf=0.8
y0=0.13 
yf=0.83
ym=0.5*(y0+yf)
;elabel=string(e0,flab)+' - '+string(e1,flab)
mtitle='mu = '+string(mu0,flab)+' keV/nT,  K = '+ $
       string(k0,fklab)+' nT^1/2 R_E '+species
read,'Do you want to reduce the temporal resolution? y(yes), n(no) => ',yon
if (yon eq 'n') then begin
   plspsds=plspsd
   hours=hour
endif else begin
   read,'enter number of time sector => ',ntsec
   ntnew=ntsec*2
   dn=fix(ntime/ntsec)
   plspsds=fltarr(ntnew,nr)   &   hours=fltarr(ntnew)
   for n=0,ntsec-1 do begin
       n0=n*dn
       n1=n0+dn
       n2=n*2
       hours(n2)=hour(n0)   
       hours(n2+1)=hour(n1)-0.01
       for i=0,nr-1 do plspsds(n2,i)=(total(plspsd(n0:n1,i)))/(n1-n0+1)
       plspsds(n2+1,0:nr-1)=plspsds(n2,0:nr-1)  
   endfor
endelse
plotPsd:
x_tick=fix(max(hour))/24
contour,plspsds,hours,Lshell,yrange=[1,8],ystyle=1,xtitle='hour',ytitle='L', $
        xrange=[min(hour),max(hour)],xstyle=1,xticks=x_tick,title=mtitle,  $
        charsize=1.5,pos=[x0,y0,xf,yf],levels=lvl,c_colors=colr,/fill
xyouts,x0,0.93,fhead,size=2.,/normal

; draw color bar and label
x1=xf+0.06
x2=x1+0.03
dy=(yf-y0)/(colrmax-colrmin+1)
for i=colrmin,colrmax do begin
    y1=y0+i*dy
    y2=y1+dy*1.02
    polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=i,/normal  
endfor
if (iunit eq 1) then bar_lab='psd [c3/(MeV3 cm3)]'
if (iunit eq 2) then bar_lab='psd [s3/(kg3 m6]'
if (ilog eq 1) then bar_lab='log '+bar_lab
xyouts,x2,y0,string(fmin,'(f5.1)'),size=1.5,/normal
xyouts,x2,yf,string(fmax,'(f5.1)'),size=1.5,/normal
xyouts,x2+0.03,ym,bar_lab,alignment=0.5,orientation=90.,size=1.5,/normal
if (iplt eq 1) then iplt_lab='equatorial psd'
if (iplt eq 2) then iplt_lab=string(fix(alt),format='(i4)')+' km altitude'
xyouts,0.5*(x0+xf),yf+0.08,iplt_lab,alignment=0.5,size=1.1,/normal

; Make gif file
if (ips eq 1) then goto,newPlot
gimg = tvrd(0,0)
;write_gif,fhead+'_'+fmiddle+'.gif',gimg,red(0:255),green(0:255),blue(0:255) 

; Make ps file
set_plot,'ps'
ips=1
device,filename='plot_psd_L_Time.ps',/inches,/color,yoffset=0.5,xoffset=0.3, $
       xsize=7.,ysize=7.*y_wsize/x_wsize
goto,plotPsd

newPlot:
read,'Do you want to continue?  (y)yes, (n)no => ',yon
if (yon eq 'y') then goto,new_plot
;write_gif,fhead+'_'+fmiddle+'.gif',/close
device,/close

end
