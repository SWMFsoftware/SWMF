;-------------------------------------------------------------------------------
;
;                              plot_fls_L_Time.pro
;
; IDL procedure reads *.fls file and plot the pitch-angle local-time averaged
; flux in a L-time plot.
;
; Created on 9 November 2005 by Mei-Ching Fok, Code 612.2, NASA GSFC.
; March 29, 2017 G. Toth made into a procedure with optional arguments.
;-------------------------------------------------------------------------------

pro plot_fls_L_Time, fhead=fhead, iplt=iplt, alt=alt, ntime=ntime, ilog=ilog, $
                     ie1=ie1, ie2=ie2, iunit=iunit, fmin=fmin, fmax=fmax, $
                     interactive=interacive

  ;; Optional arguments
  ;; fhead       - filename without extension. No default.
  ;; iplt        - plot area: 1=equator, 2=ionosphere. Default is 1.
  ;; alt         - altitude for ionosphere plot in km. Default is 100.
  ;; ntime       - number of time steps in file. Has to be more than 1. No default.
  ;; ilog        - 1=logarithmic scale, 0=linear scale. Default is 1.
  ;; ie1         - index of lowest energy bin. Default is 0.
  ;; ie2         - index of highest energy bin. Default is 11.
  ;; iunit       - flux units: 1= /keV cm2 s sr, 2= /cm2 s sr. Default is 1.
  ;; fmin        - minimum value of (log)flux. Default is 0 (linear) or fmax-3 (log).
  ;; fmax        - maximum value of (log)flux. Default is maximum found in output.
  ;; interactive - interactive mode requiring user input. Default is 0 (false).

  if not keyword_set(interactive) then interactive=0

  close,/all
  ;;write_gif,/close

  ;; get colors for color table palette
  loadct,39                     ;  rainbow with white
  tvlct, red, green, blue, /get

  ;; read file name and energy, sina, lat information 
  if n_elements(fhead) eq 0 then begin
     fhead=' '
     read,'enter the file head (e.g., 2000_225_h) => ',fhead
  endif
  openr,2,fhead+'.fls'
  readf,2,rc,ir,ip,je,ig        ; je=number of energy grid, ig=number of pa grid
  nr=ir                         ; number of radial grid
  nmlt=ip                       ; number of local-time grid
  energy=fltarr(je)    &   sina=fltarr(ig)   &   lat=fltarr(ir)   &   sinao=sina
  Ebound=fltarr(je+1)
  readf,2,energy
  ;;readf,2,Ebound
  readf,2,sina
  readf,2,lat

  if n_elements(iplt) eq 0 then $
     if interactive then read,'plot flux at  1 = equator,  2 = ionosphere => ',iplt $
     else iplt = 1
  if n_elements(alt) eq 0 then begin
     alt=100.
     if (interactive and iplt eq 2) then read,'altitude (km) => ',alt
  endif
  ri=(alt+6375.)/6375.          ; ionosphere distance in RE
  rim=ri*6.375e6                ; ionosphere distance in m
  xme=7.9e15                    ; magnetic dipole moment of the Earth

  if n_elements(ntime) eq 0 then $
     read,'number of data set in time => ',ntime

  if ntime lt 2 then begin
     print,'The file ',fhead,' has to contain at least 2 snapshots to make this plot.'
     retall
  endif

  if n_elements(ilog) eq 0 then $
     if interactive then read,'flux in?  0 = linear scale, or 1 = log scale => ',ilog $
     else ilog=1

  ;; Calculate L-values
  Lshell=lat
  for i=0,ir-1 do begin
     cosl=cos(lat(i)*!pi/180.)
     Lshell(i)=rc/cosl/cosl
  endfor

  ;; define energy labels
  for k=1,je-1 do Ebound(k)=sqrt(energy(k-1)*energy(k))
  Ebound(0)=energy(0)*energy(0)/Ebound(1)
  Ebound(je)=energy(je-1)*energy(je-1)/Ebound(je-1)
  if (Ebound(je) lt 999.5) then ilab='(i3.3)'
  if (Ebound(je) ge 999.5) then ilab='(i4.4)'
  if (Ebound(je) lt 999.5) then flab='(f5.1)'
  if (Ebound(je) ge 999.5) then flab='(f6.1)'
  
  ;; setup energy and species
  species = 'e-'
  spe=strmid(fhead,strlen(fhead)-1,1) ; species from the last character
  if (spe eq 'h') then species='H+'
  if (spe eq 'o') then species='O+'
  if (spe eq 'e') then species='e-'

  ;; setup arrays
  plsflux2=fltarr(nmlt,je)
  plsfluxk=fltarr(ntime,nr,je)    
  plsflux=fltarr(ntime,nr)      
  hour=fltarr(ntime)
  fy=fltarr(ig)      &    fynew=fy    &     dmu=fy        &    cosa=fy
  cosa(*)=cos(asin(sina(*)))
  for m=0,ig-1 do begin
     if (m eq 0) then sina0=0.
     if (m gt 0) then sina0=0.5*(sina(m)+sina(m-1))
     if (m eq (ig-1)) then sina1=1.
     if (m lt (ig-1)) then sina1=0.5*(sina(m)+sina(m+1))
     dmu(m)=sqrt(1.-sina0*sina0)-sqrt(1.-sina1*sina1)
                                ;print,  dmu(m),sina0,sina0,sina1,sina1
  endfor
  dE=fltarr(je)
  for k=0,je-1 do dE(k)=Ebound(k+1)-Ebound(k)

  ;; Read fluxes and calculate pa-mlt averaged flux
  for n=0,ntime-1 do begin
     readf,2,hour1
     hour(n)=hour1
     print,n,hour1
     for i=0,nr-1 do begin
        bi=sqrt(4.-3.*ri/Lshell(i))*xme/rim^3 ; assume dipole
        for j=0,nmlt-1 do begin
           readf,2,lat1,xmlt,ro,xmlto,bo ; bo: equatorial B in T
           sqrtbobi=sqrt(bo/bi)
           sinao(*)=sina(*)*sqrtbobi
           for k=0,je-1 do begin
              readf,2,fy
              fynew=fy
              if (iplt eq 2) then fynew=interpol(fy,sina,sinao)
              plsflux2(j,k)=0.  ; pitch-angle averaged flux
              for m=0,ig-1 do begin
                 plsflux2(j,k)=plsflux2(j,k)+fynew(m)*dmu(m)
                                ;print,plsflux2(j,k),fynew(m),dmu(m)
              endfor
           endfor
        endfor                  ; plasfluxk: pa-mlt averaged flux
        for k=0,je-1 do begin
           plsfluxk(n,i,k)=total(plsflux2(*,k))/nmlt 
        endfor
     endfor
  endfor
  close,2

new_plot:
  ;; choose which energy to be displaced
  print,'energy bin (keV): '
  for i=0,je-1 do print,i,Ebound(i:i+1),format='("  (",i2,") ",f6.1," - ",f6.1)'

  if n_elements(ie1) eq 0 or n_elements(ie2) eq 0 then $
     if interactive then read, 'lower and upper energy bins (ie1,ie2) => ',ie1,ie2 $
     else ie1=0 & ie2=11

  if n_elements(iunit) eq 0 then $
     if interactive then read,'unit? (1)flux/(keV cm2 s sr), (2)flux/(cm2 s sr) => ',iunit $
     else iunit=1
  e0=Ebound(ie1)
  e1=Ebound(ie2+1)
  fmiddle=string(round(e0),ilab)+'-'+string(round(e1),ilab)+'keV'

  ;; calculate total flux 
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
        if (iunit eq 1) then plsflux(n,i)=plsflux(n,i)/weiSum ; diff flux
     endfor
  endfor
  if (ilog eq 1) then plsflux=alog10(plsflux)

  ;; setup plot ranges
  yon=' '
  if n_elements(fmax) eq 0 then begin
     fmax=max(plsflux)    
     print,' fmax = ',fmax
     if keyword_set(interactive) then begin
        read,' Do you want to change it? (y/n) => ',yon
        if (yon eq 'y') then read,' enter new fmax => ',fmax
     endif
  endif
  if n_elements(fmin) eq 0 then begin
     fmin=0.   
     if (ilog eq 1) then fmin=fmax-5.
     if keyword_set(interactive) then begin
        read,' Do you want to change fmin? (y/n) => ',yon
        if (yon eq 'y') then read,' enter new fmin => ',fmin
     endif
  endif

  ;; limit flux to fmin
  plsflux = plsflux > fmin

  ;; setup levels
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

  ;; Setup window
  ips=0
  set_plot,'x'
  x_wsize=570
  y_wsize=446
  window,1,xpos=200,ypos=100,xsize=x_wsize,ysize=y_wsize
  !p.color=0
  !p.background=255

  ;; Plot plasma flux 
  x0=0.1
  xf=0.8
  y0=0.13 
  yf=0.83
  ym=0.5*(y0+yf)
  elabel=string(e0,flab)+' - '+string(e1,flab)
  mtitle=elabel+' keV '+species

  plsfluxs=plsflux
  hours=hour

  if interactive then begin
     read,'Do you want to reduce the temporal resolution? y(yes), n(no) => ',yon
     if (yon eq 'y') then begin
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
     endif
  endif
plotFlux:
  x_tick=fix(max(hour))/24
  contour,plsfluxs,hours,Lshell,yrange=[1,8],ystyle=1,xtitle='hour',ytitle='L', $
          xrange=[min(hour),max(hour)],xstyle=1,xticks=x_tick,title=mtitle,  $
          charsize=1.5,pos=[x0,y0,xf,yf],levels=lvl,c_colors=colr,/fill 
  xyouts,x0,0.93,strmid(fhead,0,strlen(fhead)-2),size=2.,/normal

  ;; draw color bar and label
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

  ;; Make gif file
  if (ips eq 1) then goto,newPlot
  gimg = tvrd(0,0)
  ;;write_gif,fhead+'_'+fmiddle+'.gif',gimg,red(0:255),green(0:255),blue(0:255) 

  ;; Make ps file
  set_plot,'ps'
  ips=1
  device,filename='plot_fls_L_Time.ps',/inches,/color,yoffset=0.5,xoffset=0.3, $
         xsize=7.,ysize=7.*y_wsize/x_wsize
  goto,plotFlux

newPlot:
  if interactive then begin
     read,'Do you want to continue?  (y)yes, (n)no => ',yon
     if (yon eq 'y') then goto,new_plot
  endif
  ;;write_gif,fhead+'_'+fmiddle+'.gif',/close
  device,/close

end
