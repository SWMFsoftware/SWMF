; This routine plots a file of solar wind data as prepped for the
; Weimer model.
; INPUT: weimer_input.inp
; OUTPUT: fig_sw.ps (plot)

npwind=0L
pwind='weimer_input.inp'
;
close,33
openr,33,pwind
dum = '        '
for i=1,4 do begin
  readf,33,dum
endfor
readf,33,npwind,year,doy,vavg,tshift
;
bump=fltarr(6,npwind)
;readf,33,format='(11x,f2.0,1x,f2.0,1x,f6.3,3f14.1,3f14.2,f7.2)' $
; ,bump
readf,33,bump
close,33

hrutdw = reform(bump(0,*))
hrutd = hrutdw - tshift
doy_ace = hrutd/24.0 + doy

by   = reform(bump(1,*))
bz   = reform(bump(2,*))
vx   = reform(bump(3,*))
den  = reform(bump(4,*))
tilt = reform(bump(5,*))

drange = [0,max(den)]
vrange = [min(vx),max(vx)]

tit = 'Year '+ string(format='(i4)',year)
xtit="Day of Year"
iwin=0

set_plot,'x'
  !x.style = 1
  !p.multi=0
  re = 6370.0
  xwin = 640*1.
  ywin = 640*1.
  window,xsize=xwin,ysize=ywin,iwin

  erase
      plot,doy_ace,by,xtitle=" "  $
         ,ytitle="!3B!By!N (nT)" $
         ,title=tit            $
;         ,xrange=[xmin,xmax],yrange=[zmin,zmax]            $
;         ,ytickname=["!30"] $
;         ,xtickname=["!30"] $
;         ,thick=2 $
         ,xtickname=[" "," "," "," "," "," "," "," "] $
         ,charsize=1.4,pos=[0.2,0.77,0.90,0.95] $
         ,/noerase
      oplot,[0.,400.],[0.,0.]
;
      plot,doy_ace,bz,xtitle=" "  $
         ,ytitle="!3B!Bz!N (nT)" $
         ,title=" "            $
;         ,xrange=[xmin,xmax],yrange=[zmin,zmax]            $
;         ,ytickname=["!30"] $
;         ,xtickname=["!30"] $
;         ,thick=2 $
         ,xtickname=[" "," "," "," "," "," "," "," "] $
         ,charsize=1.4,pos=[0.2,0.57,0.90,0.75] $
         ,/noerase
      oplot,[0.,400.],[0.,0.]         
;
      plot,doy_ace,vx,xtitle=" "  $
         ,ytitle="!3V!BSW!N (km s!A-1!n)" $
         ,title=" "            $
         ,yrange=vrange            $
;         ,xrange=[xmin,xmax],yrange=[zmin,zmax]            $
;         ,ytickname=["!30"] $
         ,xtickname=[" "," "," "," "," "," "," "," "]  $
         ,charsize=1.4,pos=[0.20,0.37,0.90,0.55] $
         ,/noerase

      plot,doy_ace,den,xtitle=xtit  $
         ,ytitle='n (cm!A-3!N)' $
         ,title=" "            $
         ,yrange=drange            $
;         ,ytickname=["!30"] $
         ,charsize=1.4,pos=[0.20,0.17,0.90,0.35] $
         ,/noerase

;
; -----------------------
; set plot to postscript
; -----------------------

  !p.charthick  = 4
  !p.thick      = 2
  xth           = 4
  yth           = 4  

  set_plot,'ps'
  device,file='fig_sw.ps',yoffset=1,ysize=25

;  device,file='fig.ps',bits_per_pixel=8,/color, $
;         xsize=8,ysize=5,/inches,xoffset=.25

  erase
      plot,doy_ace,by,xtitle=" "  $
         ,ytitle="!3B!By!N (nT)" $
         ,title=tit            $
;         ,xrange=[xmin,xmax],yrange=[zmin,zmax]            $
;         ,ytickname=["!30"] $
;         ,xtickname=["!30"] $
;         ,thick=2 $
         ,xtickname=[" "," "," "," "," "," "," "," "] $
         ,charsize=1.4,pos=[0.2,0.77,0.90,0.95] $
         ,/noerase
      oplot,[0.,400.],[0.,0.]
;
      plot,doy_ace,bz,xtitle=" "  $
         ,ytitle="!3B!Bz!N (nT)" $
         ,title=" "            $
;         ,xrange=[xmin,xmax],yrange=[zmin,zmax]            $
;         ,ytickname=["!30"] $
;         ,xtickname=["!30"] $
;         ,thick=2 $
         ,xtickname=[" "," "," "," "," "," "," "," "] $
         ,charsize=1.4,pos=[0.2,0.57,0.90,0.75] $
         ,/noerase
      oplot,[0.,400.],[0.,0.]         

      plot,doy_ace,vx,xtitle=" "  $
         ,ytitle="!3V!BSW!N (km s!A-1!n)" $
         ,title=" "            $
         ,yrange=vrange            $
;         ,xrange=[xmin,xmax],yrange=[zmin,zmax]            $
;         ,ytickname=["!30"] $
         ,xtickname=[" "," "," "," "," "," "," "," "]  $
         ,charsize=1.4,pos=[0.20,0.37,0.90,0.55] $
         ,/noerase

      plot,doy_ace,den,xtitle=xtit  $
         ,ytitle='n (cm!A-3!N)' $
         ,title=" "            $
         ,yrange=drange            $
;         ,ytickname=["!30"] $
         ,charsize=1.4,pos=[0.20,0.17,0.90,0.35] $
         ,/noerase

; close the postscript file 

  device,/close

; reset the display to the screen

  set_plot,'x'

  !p.charthick   = 1
  !p.thick       = 1

end

