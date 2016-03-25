pro GLSETUP2, FILE=FILE, UseBATS=UseBATS

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
  sizemap_p= mag_info.bphi_field
  sizemap_n= mag_info.btheta_field

  neqpar =mag_info.neqpar
  eqpar  =mag_info.eqpar

  xPositive = eqpar[2]
  yPositive = eqpar[3]
  xNegative = eqpar[4]
  yNegative = eqpar[5]
  XyARCenter_D = [eqpar[6],eqpar[7]]
  nPIL = (neqpar - 8)/2
  xPIL_I = eqpar[8:7+nPIL]
  yPIL_I = eqpar[8+nPIL:7+2*nPIL]

 

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
  plots,xPositive,yPositive,/data,psym=-2,color=250
  plots,xNegative,yNegative,/data,psym=-2,color=50
;plot center of the flux rope 
  plots,XyARCenter_D[0],XyARCenter_D[1],/data,psym=-2,color=150
;show PIL
  for i=0, nPIL-1 do begin
     plots,xPIL_I(i),yPIL_I(i),psym=-1,color=200
  endfor
; Showing positive and negative spots
  contour,sizemap_p,/overplot,c_color=100
  contour,abs(sizemap_n),/overplot,c_color=100
  wait, 2
;The region size is used to cover the whole area of active region in
;order to show a zoom-in image. Shorter RegionSize for near-Limb
;regions when needed.

  RegionSize=round(long(50)*nlon/360)

  
  ;Display the zoom-in image of the active region with weighted centers and PIL.
  window,3,xs=800,ys=800,xpos=400,ypos=400

  device,decomposed=0
  loadct,0
  
  sub_x1=max([round(XyARCenter_D[0])-RegionSize/2,0])
  sub_x2=min([round(XyARCenter_D[0])+RegionSize/2,nlon-1])
  sub_y1=max([round(XyARCenter_D[1])-RegionSize/2,0])
  sub_y2=min([round(XyARCenter_D[1])+RegionSize/2,nlat-1])

  contour,br_field_show[sub_x1:sub_x2,sub_y1:sub_y2],$                        
          min=-20,max=20,charsize=3,title='CME Source Region (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Degree)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1

  loadct,39
;Showing positive and negative spots
  contour,sizemap_p[sub_x1:sub_x2,sub_y1:sub_y2],/overplot,c_color=100
  contour,abs(sizemap_n[sub_x1:sub_x2,sub_y1:sub_y2]),/overplot,c_color=100
  plots,xPositive-sub_x1,yPositive-sub_y1,$
        /data,psym=-2,color=250,symsize=3,thick=3
  plots,xNegative-sub_x1,yNegative-sub_y1,$
        /data,psym=-2,color=50,symsize=3,thick=3
  plots,XyARCenter_D[0]-sub_x1,XyARCenter_D[1]-sub_y1,/data,psym=-2,color=150,symsize=3,thick=3
  for i=0, nPIL-1 do begin
     plots,xPIL_I(i)-sub_x1,yPIL_I(i)-sub_y1,psym=-1,color=200,symsize=3,thick=3
  endfor
  wait,10
  wdelete,2
  wait,10
  wdelete, 3
  exit
end
