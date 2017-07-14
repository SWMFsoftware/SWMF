
pro contour_circle, data, lons, lats, $
                    mini = mini, maxi = maxi, sym = sym, $
                    nLevels = nLevels, $
                    no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                    pos = pos, $
                    maxrange = maxrange, title = title, nolines = nolines, $
                    colorbar = colorbar, $
                    latExtra = latExtra, lonExtra=lonExtra, $
                    nocolor = nocolor, $
                    nomaxmin = nomaxmin, thick = thick, average = average, $
                    indicatemaxmin = indicatemaxmin

  if (n_elements(thick) eq 0) then thick = 2
  if (n_elements(nocolor) eq 0) then nocolor = 0
  if (n_elements(average) eq 0) then average = 0
  if (n_elements(nomaxmin) eq 0) then nomaxmin = 0

  if (n_elements(mini) eq 0) then mini = min(data)
  if (n_elements(maxi) eq 0) then maxi = max(data)

  if (n_elements(indicatemaxmin) eq 0) then indicatemaxmin=0

  if (maxi eq mini) then begin
      if (mini eq 0.0) then maxi = 1.0 $
      else $
        if (mini lt 0.0) then mini = mini*1.1 else maxi = maxi*1.1
  endif
  if (n_elements(sym) gt 0 and mini lt 0.0) then begin
      maxi = max([maxi, abs(mini)])
      mini = -maxi
  endif

  if (max(lons) lt 30.0 and max(lons) gt 10.0) then londiv = 12.0 
  if (max(lons) lt 10.0) then londiv = !pi
  if (max(lons) gt 30.0) then londiv = 180.0
  if (abs(max(lats)) lt 2.0*!pi) then lats = lats / !dtor

  if (n_elements(maxrange) eq 0) then maxrange = 90.0 - max(lats)
  if (n_elements(pos) eq 0) then begin
      ppp = 1
      space = 0.1
      pos_space, ppp, space, sizes
      get_position, ppp, space, sizes, 0, pos
  endif
  if (n_elements(nlevels) eq 0) then nlevels = 31
  if (n_elements(nolines) eq 0) then nolines = 0

  if (nlevels eq 1) then begin
     levels = [(mini+maxi)/2]
  endif else begin
     levels = findgen(nlevels)*(maxi-mini)/(nlevels-1) + mini
  endelse

  nlats = n_elements(lats)
  nlons = n_elements(lons)

  lat2d = fltarr(nlons,nlats)
  lon2d = fltarr(nlons,nlats)

  for i=0,nlats-1 do lon2d(*,i) = lons*!pi/londiv - !pi/2.0
  for i=0,nlons-1 do lat2d(i,*) = lats

  x = (90.0-lat2d)*cos(lon2d)
  y = (90.0-lat2d)*sin(lon2d)

  if (n_elements(latExtra) gt 0) then begin
      PlotExtraPoints = 1
      xExtraPoints = (90.0-latExtra)*cos(lonExtra*!pi/londiv - !pi/2.0)
      yExtraPoints = (90.0-latExtra)*sin(lonExtra*!pi/londiv - !pi/2.0)
  endif else PlotExtraPoints = 0

  loc = where(lat2d(0,*) ge 90.0-maxrange and lat2d(0,*) le 90.0, count)

  if (count gt 0) then begin

     if (not nocolor) then begin
        
        data_chop = data(*,loc)
        lmini = where(data_chop lt levels(1),cmini)
        if (cmini gt 0) then data_chop(lmini) = levels(1)
        lmaxi = where(data_chop gt levels(nLevels-2),cmaxi)
        if (cmaxi gt 0) then data_chop(lmaxi) = levels(nLevels-2)
        contour, data_chop, x(*,loc), y(*,loc), /noerase, pos = pos, $
                 levels = levels, $
                 c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
                 /cell, /follow, xstyle = 5, ystyle = 5, $
                 xrange = [-maxrange,maxrange], $
                 yrange = [-maxrange,maxrange]
     endif
     if (not nolines) then begin
        lsub = indgen((nLevels-1)/3)*3
        contour, data(*,loc), x(*,loc), y(*,loc), /noerase, pos = pos, $
                 levels = levels(lsub), $
                 /follow, xstyle = 5, ystyle = 5, $
                 xrange = [-maxrange,maxrange], $
                 yrange = [-maxrange,maxrange], thick = thick
     endif else begin
        if (indicatemaxmin) then begin
           ds = reform(data(*,loc))
           xs = reform(x(*,loc))
           ys = reform(y(*,loc))
           l = where(ds eq max(ds))
           oplot, [xs(l(0))], [ys(l(0))], psym = 4
           l = where(ds eq min(ds))
           oplot, [xs(l(0))], [ys(l(0))], psym = 2
        endif
     endelse
     plotmlt, maxrange, $
              no00 = no00, no06 = no06, no12 = no12, no18 = no18

     if (plotExtraPoints) then begin
        oplot, xExtraPoints, yExtraPoints, psym = 4, symsize = 3, thick = 5
     endif
     
     mini_tmp = min(data(*,loc))
     maxi_tmp = max(data(*,loc))

     if (average) then begin
        area = cos(lat2d(*,loc)*!dtor)
        ave = total(data(*,loc)*area)/total(area)
        mini_tmp = mean(data(*,loc))
        maxi_tmp = ave
     endif

      if (abs(mini_tmp) gt 1000.0) or (abs(maxi_tmp) gt 1000.0) or        $
         (abs(maxi_tmp) lt 0.1) then begin
        maxs = string(maxi_tmp,format="(e9.2)")
        mins = string(mini_tmp,format="(e9.2)")
      endif else begin
        maxs = string(maxi_tmp,format="(f7.2)")
        mins = string(mini_tmp,format="(f7.2)")
      endelse

      if (not nomaxmin) then begin
         p1 = pos(0)+(pos(2) - pos(0))*0.30 * (1.0 - sin(45.0*!pi/180.0))
         p2 = pos(1)+(pos(3) - pos(1))*0.30 * (1.0 - sin(45.0*!pi/180.0))
         xyouts, p1, p2, mins, /norm, charsize = 0.9, orient = -45, alignment=0.5

         p1 = pos(2)-(pos(2) - pos(0))*0.30 * (1.0 - sin(45.0*!pi/180.0))*0.95
         xyouts, p1, p2, maxs, /norm, align=0.5, charsize = 0.9, $
                 orient = 45.0
      endif
      if (n_elements(title) gt 0) then begin

          p1 = pos(0)+(pos(2) - pos(0))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.95
          p2 = pos(3)-(pos(3) - pos(1))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.95

          xyouts, p1, p2, title, $
                /norm, alignment = 0.5, charsize = 0.8, orient = 45.0

      endif

      if (n_elements(colorbar) gt 0) then begin

          p1 = pos(2)-(pos(2) - pos(0))*0.50 * (1.0 - sin(45.0*!pi/180.0))
          p2 = pos(3)-(pos(3) - pos(1))*0.50 * (1.0 - sin(45.0*!pi/180.0))

          ctpos = [p1,p2,p1+0.01,pos(3)]
          plotct, 254, ctpos, [mini,maxi], colorbar, /right

      endif

  endif

end
