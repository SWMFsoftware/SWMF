pro GLSETUP1, FILE=FILE, UseBATS=UseBATS, CMESpeed=CMESpeed

;-----------------------------------------------------------------------
; NAME:
;   GLSETUP1
;
; OUTPUTS:
;   user choise for stdout file
;
; SYSTEM REQUIREMENTS:
; Mouse with left and right button
;
;
; KEYWORDS: 
;
;   FILE = input magnetogram file (can be FITS or SWMF format).
;   UseBATS = if set, will read BATS-R-US format (2D or 3D). Default
;             will read FITS format.

;Setup the color mode and a better IDL font.
  device,decomposed=0
  !p.font=1
  PlotRadius =1.

;Read the magnetogram
  mag_info=read_magnetogram(file,PlotRadius,UseBATS)
  nlat=mag_info.nlat
  nlon=mag_info.nlon
  longitude=mag_info.longitude
  latitude=mag_info.latitude
  br_field=mag_info.br_field
  neqpar =mag_info.neqpar
  eqpar  =mag_info.eqpar
 

;Display the magnetogram and let user interactively select the CME source region. The
;procedure to select is:
; 1. Click the CME source region of positive polarity with 'left' button of mouse
; 2. Click the CME source region of negative polarity with 'right' button of mouse
;
;Note that the user can click anywhere inside the active region. However, click closer
;to the center of the positive/negative patterns is recommended.
;
;Note the solar latitude is expressed in pixel due to the non-uniform spacing. The latitude
;is uniform in sin(latitude). This will be changed in the future to degree. 

  br_field_show=br_field
  index=where(br_field lt -20)
  br_field_show[index]=-20
  index=where(br_field gt 20)
  br_field_show[index]=20

  window,2,xs=1200,ys=1200.*float(nlat)/float(nlon)*4./3.
  loadct,0
  contour,br_field_show,min=-20,max=20,charsize=3,$
          title='SWMF Input Magnetogram (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Pixel)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1
  
  
  loadct,39
  if(neqpar ge 2) then begin
     ;plot the Earth Carrington coordinate:
     if eqpar[1] gt 0 then begin
        ;eqpar[1] is the Carrington coordinate of the Earth
        ;eqpar[0] is the Carrington coordinate of the left margin
        MapLongEarth = eqpar[1]-eqpar[0]
        ;If the left margin of the map is in the next Carrington 
        ;rotation, add 360 deg:
        if MapLongEarth lt 0 then MapLongEarth = MapLongEarth +360
        ;Our plot coordinates are in pixels, not in degrees:
        PixelEarth = (MapLongEarth/360.)*nlon
        xEarthLine = findgen(nlat)*0.+PixelEarth
        yEarthLine = findgen(nlat)
        oplot,xEarthLine,yEarthLine,color=250,linestyle=5
     endif
  endif
  !MOUSE.button = 0
  while(!MOUSE.button ne 1) do begin
     cursor,xPositive,yPositive,/data,/down
     if br_field[xPositive,yPositive] lt 0 then begin
        print,'Negative Polarity! Please Select POSITIVE Polarity!'   
        !MOUSE.button=0
     endif else begin
        plots,xPositive,yPositive,/data,psym=-2,color=250 
     endelse
  endwhile
  while(!MOUSE.button ne 4) do begin
     cursor,xNegative,yNegative,/data,/down
     if br_field[xNegative,yNegative] gt 0 then begin
        print,'Positive Polarity! Please Select NEGATIVE Polarity!'   
        !MOUSE.button=0
     endif else begin
        plots,xNegative,yNegative,/data,psym=-2,color=50
     endelse
  endwhile
  print, '==='
  print, xPositive,yPositive,xNegative,yNegative,CMESpeed

  !mouse.button=0
  wait,2
  wdelete,2
  exit
end
