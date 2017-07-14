;This program can do vectors. User select factor from the slider bar.
;Select vector step from the slider bar. Tuesday 29 April 2003

pro read_thermosphere_file, filelist, nvars, nalts, nlats, nlons, vars, data, rb, cb, bl_cnt

  print, filelist(0)

  if (strpos(filelist(0), "save") gt 0) then begin
    restore, filelist(0)
    bl_cnt = 1
    return
  endif

  f = filelist(0)
  all = findfile('*'+strmid(f,6,strlen(f)-6))
;  all = findfile('*'+strmid(filelist(0),6,18))

  nfiles = n_elements(all)

  nBLKlat = fix(strcompress(cb, /REMOVE_ALL))
  nBLKlon = fix(strcompress(rb, /REMOVE_ALL))

  if (nBLKlat eq 0 and nBLKlon eq 0) then begin

    if (nfiles eq 1) then begin
      nBLKlon = 1
      nBLKlat = 1
    endif

    if (nfiles eq 2) then begin
      nBLKlon = 2
      nBLKlat = 1
    endif

    if (nfiles eq 8) then begin
      nBLKlon = 4
      nBLKlat = 2
    endif

    if (nfiles eq 16) then begin
      nBLKlon = 4
      nBLKlat = 4
    endif

    if (nfiles eq 32) then begin
      nBLKlon = 8
      nBLKlat = 4
    endif

 endif

 if (nBLKlat*nBLKlon eq 0) then begin
   bl_cnt = 0
   mes=widget_message('Please enter the block numbers! lon x lat !')
 endif else begin
  if (nBLKlat*nBLKlon gt nfiles) then begin
    bl_cnt = 0
    mess = 'There are not enough files to fill the blocks! Blocks:'+ $
           string(nBLKlat*nBLKlon)+'  files:'+string(nfiles)
    mes=widget_message(mess)
  endif else begin
    bl_cnt = 1
    for n=0,nBLKlat*nBLKlon-1 do begin
      file = all(n)
      print,'-----------------------------------------'
      print, "Reading file ", file
      openr,1,file
      done = 0
      line = ''
      while not done do begin
          readf,1,line
          if strpos(line,'NUMERICAL') gt -1 then begin
              nvars = 0L
              nlines = 0L
              readf,1, nvars
              readf,1, nalts
              readf,1, nlats
              readf,1, nlons
          endif
          if strpos(line,'VARIABLE') gt -1 then begin
              if n_elements(nvars) eq 0 then begin
                  print, 'File is in the wrong order, NUMERICAL VALUES must come'
                  print, 'before VARIABLE LIST'
                  stop
              endif else begin
                  vars = strarr(nvars)
                  for i=0,nvars-1 do begin
                      readf,1,format="(I7,a)",j,line
                      vars(i) = line
                  endfor
              endelse
              ;Store the value of vars into the pointer. 
          endif
          if strpos(line,'BEGIN') gt -1 then done = 1
      endwhile
      if (n eq 0) then begin
          nlo = (nlons-4) * nBLKlon + 4
          nla = (nlats-4) * nBLKlat + 4
          data = fltarr(nvars, nlo, nla, nalts)
          tmp = fltarr(nvars)
          format = '('+tostr(nvars)+'E11.3)'
      endif
      
;      print,''
;      print,'nalts,nlats,nlons:',nalts,nlats,nlons
;      print,'nlo,nla:',nlo,nla
;      print,'(n mod nBLKlon)*(nlons-4)',(n mod nBLKlon)*(nlons-4)
;      print,'(n/nBLKlon)*(nlats-4)=',(n/nBLKlon)*(nlats-4)
      
      for k = 0, nalts-1 do begin
          for j = 0, nlats-1 do begin
              for i = 0, nlons-1 do begin
                  readf,1,tmp, format=format
                  ii = (n mod nBLKlon)*(nlons-4) + i
		  jj = (n / nBLKlon)*(nlats-4) + j

                  ; Lets try to not put in the ghost cells...
                  if ((j ge 2 and j le nlats-3) and $
                      (i ge 2 and i le nlons-3)) then begin
                      data(*,ii,jj,k) = tmp
                  endif

                  if (j lt 2 and jj lt 2) then data(*,ii,jj,k) = tmp
                  if (j gt nlats-3 and jj gt nla-3) then data(*,ii,jj,k) = tmp

                  if (i lt 2 and ii lt 2) then data(*,ii,jj,k) = tmp
                  if (i gt nlons-3 and ii gt nlo-3) then data(*,ii,jj,k) = tmp

              endfor
          endfor	
      endfor
      close,1
    endfor
  	nlons = nlo
  	nlats = nla

;    nvars = nvars+2
;    datanew = fltarr(nvars,nlons,nlats,nalts)
;    vars = [vars, "Vn (Hor. Mag.)", "Vi (Hor. Mag.)"]
;
;    datanew(0:nvars-3, *, *, *) = data(0:nvars-3, *, *, *)
;
;    ivne = 11
;    ivnn = 12
;
;    ivie = 23
;    ivin = 24
;
;    datanew(nvars-2,*,*,*) = $
;      sqrt(datanew(ivne,*,*,*)^2+datanew(ivnn,*,*,*)^2)
;    datanew(nvars-1,*,*,*) = $
;      sqrt(datanew(ivie,*,*,*)^2+datanew(ivin,*,*,*)^2)
;
;    data = datanew

  endelse
endelse

  data = data(*,2:nlo-2,2:nla-3,*)
  nlo = nlo - 3
  nla = nla - 4
  data(*,nlo-1,*,*) = data(*,0,*,*)

end

pro plotvectors,vars,k,data, $
                lat,lon, nlats,nlons, $
                cf,vs,vi_cnt,vn_cnt,step, $
                polar, maxran
  ;Calculate the factor
  print,'cf,vs,step:',cf,vs,step
  count = 0
  vieast_index = 0
  vneast_index = 0
  if (vi_cnt eq 1) then begin
      count = n_elements(vars)
      test = "Vi (east)"
      for i=0,count-1 do begin
          var=strcompress(vars[i],/remove_all)
          tes = strcompress(test,/remove_all)
          result= STRCMP( var, tes, 8)
          if (result eq 1) then vieast_index = i
      endfor
  endif
  if (vn_cnt eq 1) then begin
      count = n_elements(vars)
      test = "Vn (east)"
      for i=0,count-1 do begin
          var=strcompress(vars[i],/remove_all)
          tes = strcompress(test,/remove_all)
          result= STRCMP( var, tes,8 )
          if (result eq 1) then vneast_index = i
      endfor
  endif

  factors = [1.0, 5.0, 10.0, 20.0, 25.0, $
             50.0, 75.0, 100.0, 150.0, 200.0]

  factor = factors(cf-1)

;  if (vi_cnt eq 1) then begin
;      factor = max(data(vieast_index,*,*,k))/cf / 2.0
;  endif else begin
;      factor = max(data(vneast_index,*,*,k))/cf / 2.0
;  endelse

  if (vi_cnt eq 1) then vindex = vieast_index  $
  else vindex = vneast_index

  for i =0,nlats-1,step do begin

      la = lat(0,i,k)
      if (90-la lt maxran) then begin

          for j =0,nlons-1,step do begin
              lo = lon(j,i,k)

              ux = data(vindex,j,i,k)/factor
              uy = data((vindex+1),j,i,k)/factor

              if (polar) then begin
                  x = (90.0 - la) * cos(lo*!pi/180.0 - !pi/2.0)
                  y = (90.0 - la) * sin(lo*!pi/180.0 - !pi/2.0)

                  ulo = ux
                  ula = uy

                  ux = - ula * cos(lo*!pi/180.0 - !pi/2.0)  $
                       - ulo * sin(lo*!pi/180.0 - !pi/2.0)
                  uy = - ula * sin(lo*!pi/180.0 - !pi/2.0) $
                       + ulo * cos(lo*!pi/180.0 - !pi/2.0)

              endif else begin
                  x = lo
                  y = la
              endelse

              ;ux is the eastward welocity (neutral or ion)
              ;uy is the northward velocity (neutral or ion)
              oplot,[x],[y],psym = 4, color = 0
              oplot,[x,x+ux],[y,y+uy], color = 0

          endfor

      endif

  endfor

  if (polar) then begin

      x  =  maxran
      y  =  maxran
      uy = -10.0
      ux =   0.0

      plots,[x],[y],psym = 4
      plots,[x,x+ux],[y,y+uy]

      str = string(factor*10.0, format = '(f6.1)') + ' m/s'
      xyouts, x-1.0, y+uy/2.0, str, alignment = 1.0

  endif else begin

      x  =  360.0
      y  =   95.0
      uy =    0.0
      ux =  -10.0

      plots,[x],[y],psym = 4
      plots,[x,x+ux],[y,y+uy]

      str = string(factor*10.0, format = '(f6.1)') + ' m/s'
      xyouts, x, y+5.0, str, alignment = 1.0

  endelse

end


pro thermosphere_vector2_event,event

widget_control,event.top,get_uvalue=ptr,/no_copy
widget_control,event.id,get_uvalue=whichevent

case whichevent of 
	'ROWBLOCK':begin
		widget_control,(*ptr).rb_txt,get_value=txt
		(*ptr).rb = txt[0]
		print,(*ptr).rb
		end
	'COLBLOCK':begin
		widget_control,(*ptr).cb_txt,get_value=txt
		(*ptr).cb = txt[0]
		print,(*ptr).cb
		end
	'FILENAME':begin
		widget_control,(*ptr).file_txt,get_value=txt
		(*ptr).filename=txt[0]
		print,(*ptr).filename
	end
	'SELVARS':begin
		sel=event.index
		(*ptr).sel=sel
		print,'sel',sel
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt
	end
	'SEARCH':begin
		(*ptr).cursor_cnt=0
		widget_control,/hourglass
		;Get the file
		widget_control,(*ptr).file_txt,get_value=txt
		(*ptr).filename=txt[0]

		widget_control,(*ptr).cb_txt,get_value=txt
		(*ptr).cb = txt[0]
		widget_control,(*ptr).rb_txt,get_value=txt
		(*ptr).rb = txt[0]
		rb = (*ptr).rb
		cb = (*ptr).cb
		
		print,'rb & cb:',(*ptr).rb,(*ptr).cb
		file=(*ptr).filename
		if file eq '' then begin
			mes=widget_message('Please enter the file to search!')
		endif else begin
			widget_control,/hourglass
			print,'file is:',file
			filelist = findfile(file)
			nfiles = n_elements(filelist)
			(*ptr).nfiles=nfiles
			print,'filelist is:',filelist
			print,'nfiles:',nfiles
		  if filelist[0] eq '' or nfiles eq 0 then begin
			me=widget_message('There is no such a file!')
		  endif else begin
                      print, "going into read_thermo..."
                      read_thermosphere_file, filelist, nvars, nalts, nlats, nlons,vars,data,rb,cb,bl_cnt
			print,'bl_cnt:',bl_cnt
			(*ptr).bl_cnt = bl_cnt
		    if (*ptr).bl_cnt eq 1 then begin
                      ;Store the array 'data' to the '*(*ptr).data'
                      *(*ptr).vars=vars
                      *(*ptr).data=data

                      ;Store the nvars,nalts,nlats,nlons to the pointer.
                      (*ptr).nvars=nvars
                      (*ptr).nalts=nalts
                      (*ptr).nlats=nlats
                      (*ptr).nlons=nlons
                      print,'nvars,nalts,nlats,nlons:'
                      print,nvars,nalts,nlats,nlons
		    endif
		  endelse
		endelse	
	      if (*ptr).bl_cnt eq 1 then begin
		;Set the VARS to the vars_list
		widget_control,(*ptr).vars_list,set_value=*(*ptr).vars
		
		;Reset the maximum number of the sliders
		cnt1=(*ptr).lat_lon_cnt
		cnt2=(*ptr).alt_lon_cnt
		cnt3=(*ptr).alt_lat_cnt
		if cnt1 eq 1 then begin
		  widget_control,(*ptr).sli,set_slider_max=nalts-1
		endif
		if cnt2 eq 1 then begin
		  widget_control,(*ptr).sli,set_slider_max=nlats-1
		endif
		if cnt3 eq 1 then begin
		  widget_control,(*ptr).sli,set_slider_max=nlons-1
		endif
	      endif
	end
	'LAT_LON':begin
		(*ptr).cursor_cnt=0
		widget_control,(*ptr).sli,sensitive=1,set_slider_max=(*ptr).nalts-1
		(*ptr).lat_lon_cnt=1
		(*ptr).alt_lon_cnt=0
		(*ptr).alt_lat_cnt=0
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt

		widget_control,(*ptr).sli,get_value=max
	        data=*(*ptr).data			;Get the array data which be plotted. 
		mv = data(2,0,0,max)/1000.0
	        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')

		end
	'ALT_LON':begin
		(*ptr).cursor_cnt=0
		widget_control,(*ptr).sli,sensitive=1,set_slider_max=(*ptr).nlats-1
		(*ptr).lat_lon_cnt=0
		(*ptr).alt_lon_cnt=1
		(*ptr).alt_lat_cnt=0
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt

		widget_control,(*ptr).sli,get_value=max
	        data=*(*ptr).data			;Get the array data which be plotted. 
		mv = data(1,0,max,0)*180.0/!pi
	        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')

		end
	'ALT_LAT':begin
		(*ptr).cursor_cnt=0
		widget_control,(*ptr).sli,sensitive=1,set_slider_max=(*ptr).nlons-1
		(*ptr).lat_lon_cnt=0
		(*ptr).alt_lon_cnt=0
		(*ptr).alt_lat_cnt=1
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt

		widget_control,(*ptr).sli,get_value=max
	        data=*(*ptr).data			;Get the array data which be plotted. 
		mv = data(0,max,0,0)*180.0/!pi
	        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')

		end
	'SLI':begin
		(*ptr).cursor_cnt=0
		widget_control,(*ptr).sli,get_value=max
		;Get value from user regarding the button they chose.
		cnt1=(*ptr).lat_lon_cnt
		cnt2=(*ptr).alt_lon_cnt
		cnt3=(*ptr).alt_lat_cnt
	    	data=*(*ptr).data			;Get the array data which be plotted. 
		if cnt1 eq 1 then begin
			(*ptr).lat_lon_max=max
			mv = data(2,0,0,max)/1000.0
		endif
		if cnt2 eq 1 then begin
			(*ptr).alt_lon_max=max
			mv = data(1,0,max,0)*180.0/!pi
		endif
		if cnt3 eq 1 then begin
			(*ptr).alt_lat_max=max
			mv = data(0,max,0,0)*180.0/!pi
		endif
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt
	        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')
	end
	'PLOT':begin
	   ;Get cursor_x & cursor_y from widget_text.
	   widget_control,(*ptr).curx_txt,get_value=txt
	   strx=txt[0]
	   cursor_x = double(strx)
	   widget_control,(*ptr).cury_txt,get_value=txt
	   stry=txt[0]
	   cursor_y = double(stry)
	   print,'cursor_x,y:',cursor_x,cursor_y
	    
	   ;Set the value of step. 
	   if ((*ptr).vyes_cnt eq 1) then begin
		step = (*ptr).vs
		print,'step is:',step
	   endif
	 
	   ;Clean up the display window
	   plotdumb

	   ;Initialize values
	   vars=*(*ptr).vars  	   ;Get the variables' list
	   sel=(*ptr).sel	   ;Get the selected variable.
	   nfiles=(*ptr).nfiles	   ;Get the files users searched for 
	   cnt1=(*ptr).lat_lon_cnt ;Get the counter of selecting lat_lon button
           cnt2=(*ptr).alt_lon_cnt ;Get the counter of selecting alt_lon button
	   cnt3=(*ptr).alt_lat_cnt ;Get the counter of selecting alt_lat button 
	   yes=(*ptr).yes_cnt	   ;Get the counter of selecting yes (Ghost Cells) button 
	   no=(*ptr).no_cnt	   ;Get the counter of selecting no (Ghost Cells) button
	   yeslog=(*ptr).yes_logcnt;Get the counter of selecting yes Plot Log button
	   nolog=(*ptr).no_logcnt  ;Get the counter of selecting no Plot Log button
	   nalts=(*ptr).nalts	   ;Get the number of nalts
	   nlats=(*ptr).nlats	   ;Get the number of nlats
	   nlons=(*ptr).nlons	   ;Get the number of nlons

           polar = (*ptr).polar

	   showgridyes = (*ptr).gridyes_cnt
	   plotvectoryes = (*ptr).vyes_cnt
	   vi_cnt = (*ptr).vi_cnt
	   vn_cnt = (*ptr).vn_cnt
	   cf = (*ptr).cf
	   vs = (*ptr).vs
	   ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
	   cnt=1
	   if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
		cnt=0
	   endif
	   
	   ;Check if there are any variables in the text_list. 
	   ;If there is not, then initialize the counter of selecting 
	   ;plot button (plot_cnt) to zero and pop out a warning message. 
	   ;If there are variables, then start to do the plot.  
	   if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
		(*ptr).plot_cnt=0
		me=widget_message('Please do search, select a variable and a setting first!')
	   endif else begin
	    (*ptr).plot_cnt=1			;Set the plot counter to 1.
	    ;Initialize the data,alt,lat,lon,k,xrange,yrange
	    data=*(*ptr).data			;Get the array data which be plotted. 
	    alt=reform(data(2,*,*,*))/1000.0	;Calculate alt,lat,lon and k.
	    lat=reform(data(1,*,*,*))*180.0/!pi
	    lon=reform(data(0,*,*,*))*180.0/!pi
	    k=n_elements(alt(0,0,*))-1
	    xrange=[0,0]			;Initialize xrange and yrange.
	    yrange=[0,0]

	    ;According to the button user selected, get the values for plotting
	    if cnt1 eq 1 then begin
		k=(*ptr).lat_lon_max

                if (polar) then MinLat = 40.0 else MinLat = -1000.0
                mr = 90.0 - MinLat

                loc = where(lat(0,*,0) ge MinLat)

		datatoplot=reform(data(sel,*,loc,k))
		maxi=max(datatoplot)
		mini=min(datatoplot)

                if (polar) then begin

                    if (strpos((*ptr).filename,"dat") gt 0) then begin
                        hour = float(strmid((*ptr).filename,14, 2))
                        minu = float(strmid((*ptr).filename,16, 2))
                        seco = float(strmid((*ptr).filename,18, 2))
                    endif else begin
                        hour = float(strmid((*ptr).filename, 8, 2))
                        minu = float(strmid((*ptr).filename,10, 2))
                        seco = float(strmid((*ptr).filename,12, 2))
                    endelse

                    ut = hour + minu/60.0 + seco/3600.0
                    utrot = ut * 15.0 ; put into degrees

                    x = reform( (90.0 - lat(*,loc,k)) * $
                                cos((lon(*,loc,k)+utrot)*!pi/180.0 - !pi/2.0))
                    y = reform( (90.0 - lat(*,loc,k)) * $
                                sin((lon(*,loc,k)+utrot)*!pi/180.0 - !pi/2.0))
                    xrange = [-mr,mr]
                    yrange = [-mr,mr]
                    xtitle=' '
                    ytitle=' '

                endif else begin

                    utrot = 0.0

                    x=reform(lon(*,loc,k))
                    y=reform(lat(*,loc,k))
                    if yes eq 1 then begin
                        xrange=mm(lon)
                        yrange=mm(lat)
                    endif
                    if yes eq 0 then begin
                        xrange=[0,360]
                        yrange=[-90,90]
                    endif
                    xtitle='Longitude (deg)'
                    ytitle='Latitude (deg)'
                endelse
                ygrids = n_elements(loc)
                xgrids = nlons
                location = string(alt(0,0,k),format='(f5.1)')+' km Altitude'
	     endif
	     if cnt2 eq 1 then begin
		j=(*ptr).alt_lon_max
		datatoplot=reform(data(sel,*,j,*))
		maxi=max(datatoplot)
		mini=min(datatoplot)
		x=reform(lon(*,j,*))
		y=reform(alt(*,j,*))
                location = string(lat(0,j,0),format='(f5.1)')+' deg Latitude'
		xtitle='Longitude (deg)'
		ytitle='Altitude (km)'
		if yes eq 1 then begin
		   xrange=mm(lon)
		   yrange=mm(alt)
		endif
		if yes eq 0 then begin
		   backup_xrange=mm(lon)
		   backup_yrange=mm(alt)
		   default_xrange=[0,360]
		   default_yrange=mm(alt)
		  ;If out of range then use 'mm' to set xrange and yrange values.
		  ;Else use default values.
		  if (backup_xrange[0] lt default_xrange[0]) $
			and (backup_xrange[1] gt default_xrange[1]) then begin
		  	xrange=mm(lon)
			yrange=mm(alt)
		  endif else begin
			xrange=[0,360]
		   	yrange=mm(alt)
		  endelse
		endif
		ygrids=nalts
		xgrids=nlons
	     endif
	     if cnt3 eq 1 then begin
		i=(*ptr).alt_lat_max
		datatoplot=reform(data(sel,i,*,*))
		maxi=max(datatoplot)
		mini=min(datatoplot)
		x=reform(lat(i,*,*))
		y=reform(alt(i,*,*))
                location = string(lon(i,0,0),format='(f5.1)')+' deg Longitude'
		xtitle='Latitude (deg)'
		ytitle='Altitude (km)'
		if yes eq 1 then begin
		   xrange=mm(lat)
		   yrange=mm(alt)
		endif
		if yes eq 0 then begin
		  backup_xrange=mm(lat)
		  backup_yrange=mm(alt)
		  default_xrange=[-90,90]
		  default_yrange=mm(alt)
		  ;If out of range then use 'mm' to set xrange and yrange values.
		  ;Else use default values.
		  if (backup_xrange[0] lt default_xrange[0]) $
			and (backup_xrange[1] gt default_xrange[1]) then begin
		  	xrange=mm(lat)
			yrange=mm(alt)
		  endif else begin
			xrange=[-90,90]
		   	yrange=mm(alt)
		  endelse
		endif
		ygrids=nalts
		xgrids=nlats
	     endif

	   ;Calculate the xld, yld according to the cursor position user set.
	   ;Calculate and get the array will be plotted.  
	   xld = x(*,0)
	   yld = y(0,*)
	   dist_x=abs(xld-cursor_x)
	   dist_y=abs(yld-cursor_y)
	   locx=where(dist_x eq min(dist_x))
	   locy=where(dist_y eq min(dist_y))
	   datald=reform(data(sel,*,locx,locy))

	   widget_control,(*ptr).min_txt,get_value=smini
	   widget_control,(*ptr).max_txt,get_value=smaxi

	   if (float(smini) ne 0 or float(smaxi) ne 0) then begin
             mini = float(smini)
             maxi = float(smaxi)
	     mini = mini(0)
	     maxi = maxi(0)
	   endif 

	   print,'mini,maxi out:',mini,maxi
	   if mini eq maxi then maxi=mini*1.01+1.0
	   levels = findgen(31)/30.0*(maxi-mini) + mini
           loc = where(datatoplot lt levels(1), count)
           if (count gt 0) then datatoplot(loc) = levels(1)

	   ; Check if user wants to write the result to a file
	   ;If user wanted then setdevice. 
	   if (*ptr).yeswrite_cnt eq 1 then begin
		widget_control,(*ptr).writefile_txt,get_value=txt
		file_name=strtrim(string(txt[0]),2)
		;Check if user has already entered a file name. 
		if (file_name eq '') or (n_elements(file_name) eq 0) then begin
			mes=widget_message('Please enter the PS file name!')
		endif else begin
			setdevice,file_name,'l',4,0.95
		endelse
	   endif

           vars = *(*ptr).vars
           variable = strcompress(vars(sel),/remove)

	   ;Set colors. 
	   ;readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
           makect,'mid'
	   clevels = findgen(31)/30 * 253.0 + 1.0

           if (polar) then begin
               xstyle = 5 
               ystyle = 5 
           endif else begin
               xstyle = 1
               ystyle = 1
           endelse

           if (not polar and cnt1) then ppp = 2 else ppp = 1

           space = 0.075
           pos_space, ppp, space, sizes, ny = 1
           get_position, ppp, space, sizes, 0, pos
           if (not polar) then begin
               if (cnt1) then begin
                   r = pos(2) - pos(0)
                   pos(2) = pos(0) + r*2.0
               endif else begin
                   get_position, ppp, space, sizes, 0, pos, /rect
                   pos(0) = pos(0) + space
                   pos(2) = pos(2) - space
               endelse
           endif

	   ;If user DOES NOT do the Ghost Cells Contour. 
           if yes eq 0 then begin

	       ;If user DO NOT do plot log. 
               if nolog eq 1 then begin

                  ;maxi=max(datatoplot)
                  ;mini=min(datatoplot)
                   if mini eq maxi then maxi=mini*1.01+1.0
                   levels = findgen(31)/30.0*(maxi-mini) + mini

                   contour,datatoplot(*,*),x(*,*),y(*,*),POS=pos,$
                     levels=levels,xstyle=xstyle,ystyle=ystyle,$
                     xrange=xrange,yrange=yrange,$
                     title=variable+' at '+location,c_colors=clevels,$
                     xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE
                
                   if (polar) then plotmlt, mr, /no06, /no12

                   if (cnt1 eq 1) then begin
                       if (plotvectoryes eq 1) then begin
                           k=(*ptr).lat_lon_max
                           plotvectors,vars,k,data,lat,lon+utrot, $
                             nlats,nlons,cf,vs,vi_cnt,vn_cnt,step, polar, mr
                       endif 
                   endif
	      

                  ;Draw grid.
                   if (showgridyes eq 1) then begin
                       for i=0,ygrids-1 do begin
                           oplot,x(*,i),y(*,i)
                       endfor
                       for j=0,xgrids-1 do begin
                           oplot,x(j,*),y(j,*)
                       endfor
                   endif
		;If user set cursor position to plot, then do plot of datald.
		;Else clean up the cursor text fields on the interface. 
                   if (*ptr).cursor_cnt eq 1 then begin
                       if cnt1 eq 1 then begin
                           print,'do plot for Lat&Lon....'
                           plot,datald,alt(0,0,*),xtitle='datald',ytitle='Alt. (deg)'
                       endif
                       if cnt2 eq 1 then begin
                           print,'do plot for Alt&Lon........'
                           plot,datald,lat(0,*,0),ystyle=1,xtitle='datald',ytitle='Lat. (deg)'
                       endif
                       if cnt3 eq 1 then begin
                           print,'do plot for Alt&Lat....'
                           plot,datald,lon(*,0,0),xtitle='datald',ytitle='Long. (deg)'
                       endif
                   endif else begin	
                       txt=''
                       widget_control,(*ptr).curx_txt,set_value=txt
                       widget_control,(*ptr).cury_txt,set_value=txt
                   endelse
               endif            ;End if nolog eq 0.

	     ;If user DO plot log. 
               if yeslog eq 1 then begin
                   gtz = where(datatoplot gt 0,cnt)
                   if (cnt gt 0) then mini = min(datatoplot(gtz)) else mini = 1.0
                   nzsubs=where(datatoplot le 0, cnt)
                   if cnt gt 0 then datatoplot(nzsubs)=mini
                   datatoplot=ALOG10(datatoplot)
                   maxi=max(datatoplot)
                   mini=min(datatoplot)
                   print,'maxi,mini:',maxi,mini
                   if mini eq maxi then maxi=mini*1.01+1.0
                   levels = findgen(31)/30.0*(maxi-mini) + mini

                   contour,datatoplot(*,*),x(*,*),y(*,*),POS=pos,$
                     levels=levels,xstyle=xstyle,ystyle=ystyle,$
                     xrange=xrange,yrange = yrange,$
                     title=variable+' at '+location,c_colors=clevels,$
                     xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE
                   
                   if (polar) then plotmlt, mr, /no06, /no12

                   if cnt1 eq 1 then begin
                       if (plotvectoryes eq 1) then begin
                           plotvectors,vars,k,data,lat,lon+utrot,nlats,nlons, $
                                       cf,vs,vi_cnt,vn_cnt,step, polar, mr
                       endif
                   endif

              ;Draw grid.
                   if (showgridyes eq 1) then begin
                       for i=0,ygrids-1 do begin
                           oplot,x(*,i),y(*,i)
                       endfor
                       for j=0,xgrids-1 do begin
                           oplot,x(j,*),y(j,*)
                       endfor
                   endif

		;If user set cursor position to plot, then do plot of datald.
		;Else clean up the cursor text fields on the interface. 
                   if (*ptr).cursor_cnt eq 1 then begin
                       if cnt1 eq 1 then begin
                           plot,datald,alt,xtitle='datald',ytitle='Alt. (deg)'
                       endif
                       if cnt2 eq 1 then begin
                           plot,datald,lat,xtitle='datald',ytitle='Lat. (deg)'
                       endif
                       if cnt3 eq 1 then begin
                           plot,datald,lon,xtitle='datald',ytitle='Long. (deg)'
                       endif
                   endif else begin	
                       txt=''
                       widget_control,(*ptr).curx_txt,set_value=txt
                       widget_control,(*ptr).cury_txt,set_value=txt
                   endelse		
               endif            ;End if yeslog eq 1
           endif                ;ENd if yes eq 0

	   ;If user DOES the Ghost Cells Contour.
           if yes eq 1 then begin
               x1=min(lon)
               x2=max(lon)
               y1=min(lat)
               y2=min(lat)
               print,'x1,x2,y1,y2:',x1,x2,y1,y2
	     ;If user does not want to do the Plot Log.
               if nolog eq 1 then begin
                   maxi=max(datatoplot)
                   mini=min(datatoplot)
                   print,'maxi,mini 2:',maxi,mini
                   if mini eq maxi then maxi=mini+1
                   levels = findgen(31)/30.0*(maxi-mini) + mini

                   contour,datatoplot(*,*),x(*,*),y(*,*),POS=pos,$
                     levels=levels,xstyle=1,ystyle=1,$
                     xrange=xrange,yrange=yrange,$
                     title='Contour Plot Thermosphere',c_colors=clevels,$
                     xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE
	 
                   if cnt1 eq 1 then begin	
                       if (plotvectoryes eq 1) then begin
                           plotvectors,vars,k,data,lat,lon,nlats,nlons,cf,vs,vi_cnt,vn_cnt,step
                       endif
                   endif
                   
                                ;Draw grid.
                   if (showgridyes eq 1) then begin
                       for i=0,ygrids-1 do begin
                           oplot,mm(x),[y(0,i),y(0,i)]
                       endfor
                       for j=0,xgrids-1 do begin
                           oplot,[x(j,0),x(j,0)],mm(y)
                       endfor
                   endif

		;If user set cursor position to plot, then do plot of datald.
                ;Else clean up the cursor text fields on the interface. 
                   if (*ptr).cursor_cnt eq 1 then begin
                       if cnt1 eq 1 then begin
                           plot,datald,alt,xtitle='datald',ytitle='Alt. (deg)'
                       endif
                       if cnt2 eq 1 then begin
                           plot,datald,lat,xtitle='datald',ytitle='Lat. (deg)'
                       endif
                       if cnt3 eq 1 then begin
                           plot,datald,lon,xtitle='datald',ytitle='Long. (deg)'
                       endif
                   endif else begin	
                       txt=''
                       widget_control,(*ptr).curx_txt,set_value=txt
                       widget_control,(*ptr).cury_txt,set_value=txt
                   endelse
               endif            ;End of if nolog eq 1

	     ;If user does want to do the Plot Log. 
               if yeslog eq 1 then begin
                   nzsubs=where(datatoplot gt 0, cnt)
                   datatoplot(nzsubs)=datatoplot(nzsubs)
                   datatoplot=ALOG10(datatoplot)
                   maxi=max(datatoplot)
                   mini=min(datatoplot)
                   if mini eq maxi then maxi=mini+1
                   levels = findgen(31)/30.0*(maxi-mini) + mini	
                   
                   contour,datatoplot(*,*),x(*,*),y(*,*),POS=pos,$
                     levels=levels,xstyle=1,ystyle=1,$
                     xrange=xrange,yrange=yrange,$
                     title='Contour Plot Thermosphere',c_colors=clevels,$
                     xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE

                   if cnt1 eq 1 then begin
                       if (plotvectoryes eq 1) then begin
                           plotvectors,vars,k,data,lat,lon,nlats,nlons,cf,vs,vi_cnt,vn_cnt,step
                       endif
                   endif
                   
		;Draw grid.
                   if (showgridyes eq 1) then begin
                       for i=0,ygrids-1 do begin
                           oplot,mm(x),[y(0,i),y(0,i)]
                       endfor
                       for j=0,xgrids-1 do begin
                           oplot,[x(j,0),x(j,0)],mm(y)
                       endfor
                   endif
                   print,'Show Ghost Cells and do plot log'
                   print,'xrange,yrange:',xrange,yrange
                   print,'xgrids,ygrids:',xgrids,ygrids
		
		;If user set cursor position to plot, then do plot of datald.
		;Else clean up the cursor text fields on the interface. 
                   if (*ptr).cursor_cnt eq 1 then begin
                       if cnt1 eq 1 then begin
                           plot,datald,alt,xtitle='datald',ytitle='Alt. (deg)'
                       endif
                       if cnt2 eq 1 then begin
                           plot,datald,lat,xtitle='datald',ytitle='Lat. (deg)'
                       endif
                       if cnt3 eq 1 then begin
                           plot,datald,lon,xtitle='datald',ytitle='Long. (deg)'
                       endif
                   endif else begin	
                       txt=''
                       widget_control,(*ptr).curx_txt,set_value=txt
                       widget_control,(*ptr).cury_txt,set_value=txt
                   endelse
               endif            ;End of if yeslog eq 1
           endif                ;End of if yes eq 1

	   ;Draw color bar.

;	   pos = [0.82,0.05,0.87,0.96]
	   pos(0) = pos(2)+0.01
           pos(2) = pos(0)+0.03
	   maxmin = mm(levels)
	   title = ''
	   plotct,254,pos,maxmin,title,/right,color=color

	   ;If user write the result to a file, then closedevice right now. 
	   if (*ptr).yeswrite_cnt eq 1 then begin
		closedevice
		mes=widget_message('Done with writing into file!')
	   endif
	   
	   ;Show min and max number to the interface.
	   min=strtrim(string(mini),2)
	   max=strtrim(string(maxi),2)
	   widget_control,(*ptr).min_txt,set_value=min
	   widget_control,(*ptr).max_txt,set_value=max
	 endelse
	   ;Reset set cursor counter to zero.
	   (*ptr).cursor_cnt=0
	end
	'MINI':
	'MAXI':
	'GRIDYES':begin
			(*ptr).gridyes_cnt = 1
		end
	'GRIDNO':begin
			(*ptr).gridyes_cnt = 0
		end
	'VYES':begin
			widget_control,(*ptr).vector_buts,sensitive = 1
;			widget_control,(*ptr).vector_buts,/set_button
			widget_control,(*ptr).vi_but,sensitive = 1
;			widget_control,(*ptr).cf_sli,/set_button
			widget_control,(*ptr).cf_sli,sensitive=1
;			widget_control,(*ptr).vs_sli,/set_button
			widget_control,(*ptr).vs_sli,sensitive=1
			(*ptr).vyes_cnt = 1
			(*ptr).vs_cnt=1
			(*ptr).cf_cnt=1
			(*ptr).vn_cnt=1
			(*ptr).vi_cnt=0
		end
	'VNO':begin
			widget_control,(*ptr).vector_buts,sensitive = 0
;			widget_control,(*ptr).vector_buts,set_button=0
;			widget_control,(*ptr).cf_sli,/set_button
			widget_control,(*ptr).cf_sli,sensitive=0
;			widget_control,(*ptr).vs_sli,/set_button
			widget_control,(*ptr).vs_sli,sensitive=0	
			widget_control,(*ptr).vs_sli,get_value=vs
			widget_control,(*ptr).cf_sli,get_value=cf
			(*ptr).vyes_cnt = 0
			(*ptr).vs_cnt=0
			(*ptr).cf_cnt=0
			(*ptr).vs_cnt=0
			(*ptr).cf_cnt=0
			(*ptr).vs=vs
			(*ptr).cf=cf
		end
	'VI':begin
			(*ptr).vi_cnt=1
			(*ptr).vn_cnt=0
		end
	'VN':begin
			(*ptr).vi_cnt=0
			(*ptr).vn_cnt=1
		end
	'CF':begin
			widget_control,(*ptr).cf_sli,get_value=cf
			(*ptr).cf=cf
		end
	'VS':begin
			widget_control,(*ptr).vs_sli,get_value=vs
			(*ptr).vs=vs
		end
	'YES':begin
		(*ptr).cursor_cnt=0
		(*ptr).yes_cnt=1
		(*ptr).no_cnt=0
		widget_control,(*ptr).no_but,set_button=0
		widget_control,(*ptr).yes_but,sensitive=1
		widget_control,(*ptr).no_but,sensitive=1
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt
	end
	'NO':begin
		(*ptr).cursor_cnt=0
		(*ptr).no_cnt=1
		(*ptr).yes_cnt=0
		widget_control,(*ptr).yes_but,set_button=0
		widget_control,(*ptr).yes_but,sensitive=1
		widget_control,(*ptr).no_but,sensitive=1
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt
            end

        'POLAR' : begin
                   print, "Polar"
                   (*ptr).polar = 1
                 end

        'CART'  : begin
                   print, "Cartesean"
                   (*ptr).polar = 0
                 end

	'LOGYES':begin
		print,'LOG YES.'
		(*ptr).cursor_cnt=0
		(*ptr).yes_logcnt=1
		(*ptr).no_logcnt=0
		widget_control,(*ptr).no_logbut,set_button=0
		widget_control,(*ptr).no_logbut,sensitive=1
		widget_control,(*ptr).yes_logbut,sensitive=1
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt
	end
	'LOGNO':begin
		(*ptr).cursor_cnt=0
		(*ptr).yes_logcnt=0
		(*ptr).no_logcnt=1
		widget_control,(*ptr).yes_logbut,set_button=0
		widget_control,(*ptr).yes_logbut,sensitive=1
		widget_control,(*ptr).no_logbut,sensitive=1
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		widget_control,(*ptr).min_txt,set_value=txt
		widget_control,(*ptr).max_txt,set_value=txt
	end
	'SETCURSOR':begin
		(*ptr).cursor_cnt=1
		plot_cnt=(*ptr).plot_cnt
		;Reset text widgets
		txt=''
		widget_control,(*ptr).curx_txt,set_value=txt
		widget_control,(*ptr).cury_txt,set_value=txt
		if plot_cnt eq 0 then begin
		  mes=widget_message("Please do plot first!")
		endif else begin
		  cursor,x,y,/wait,/DATA
		  (*ptr).cursor_x=x
		  (*ptr).cursor_y=y
		  print,'x,y:',x,y
		  x=strtrim(string(x),2)
		  y=strtrim(string(y),2)
		  widget_control,(*ptr).curx_txt,set_value=x
		  widget_control,(*ptr).cury_txt,set_value=y
		  cnt1=(*ptr).lat_lon_cnt
           	  cnt2=(*ptr).alt_lon_cnt
	   	  cnt3=(*ptr).alt_lat_cnt
		  datatoplot=*(*ptr).data
		  help,datatoplot
		  print,'cnt1,cnt2,cnt3:',cnt1,cnt2,cnt3
		endelse
	end
	'CURX':
	'CURY':
	'YESWRITE':begin
		(*ptr).cursor_cnt=0
		(*ptr).yeswrite_cnt=1
		(*ptr).nowrite_cnt=0
		widget_control,(*ptr).nowrite_but,set_button=0
		widget_control,(*ptr).yeswrite_but,sensitive=1
		widget_control,(*ptr).nowrite_but,sensitive=1
		widget_control,(*ptr).writefile_txt,sensitive=1
		
	end
	'NOWRITE':begin
		(*ptr).cursor_cnt=0
		(*ptr).yeswrite_cnt=0
		(*ptr).nowrite_cnt=1
		widget_control,(*ptr).yeswrite_but,set_button=0
		widget_control,(*ptr).yeswrite_but,sensitive=1
		widget_control,(*ptr).nowrite_but,sensitive=1
		widget_control,(*ptr).writefile_txt,sensitive=0
	end
	'WRITEFILE':
	'CAN':begin
		widget_control,event.top,/destroy
	end
endcase
if (whichevent ne 'CAN') and (whichevent ne 'OK') then begin
	widget_control,event.top,set_uvalue=ptr,/no_copy
endif
return
NOTDONE:
	widget_control,(*ptr).list,scr_xsize=120,set_value=files
return
end

pro thermosphere_vector2

title='Plot Thermoshpere'
mainbase=widget_base(title=title,/column)
base = widget_base(mainbase,/column)
txt=''

block_base=widget_base(base,/column,/frame)
blockn_base=widget_base(block_base,/row)
lab = widget_label(blockn_base,value='Enter block number:',/align_left)
rb_txt=widget_text(blockn_base,value=txt,uvalue='ROWBLOCK',xsize=6,/editable,/align_left)
lab = widget_label(blockn_base,value='X',/align_left)
cb_txt=widget_text(blockn_base,value=txt,uvalue='COLBLOCK',xsize=6,/editable,/align_left)

dirlab=widget_base(base,/column,/frame)
lab=widget_label(dirlab,value='Enter file name (can use *) [p0000b0001i0200.dat]:',/align_left)
filename_base=widget_base(dirlab,/row)
filelist = findfile("-t *.3D* *.dat *.2D")
txt=filelist(0)
file_txt=widget_text(filename_base,value=txt,uvalue='FILENAME',xsize=46,/editable,/align_left)
sear_but=widget_button(filename_base,value='Search',uvalue='SEARCH',xsize=100)

lists_base=widget_base(base,/row,/frame)
vars_base=widget_base(lists_base,/column)
sel=-1
lab=widget_label(vars_base,value='Variables',/align_center)
vars=''
vars_list=widget_list(vars_base,value=vars,uvalue='SELVARS',xsize=20,$
	ysize=40,/align_right)
nalts=1
nlats=1
nlons=1
butt_base=widget_base(lists_base,/column)
lab=widget_label(butt_base,value='Select Settings',/align_center)
but_base=widget_base(butt_base,/column,/align_center,/frame)
buts_base=widget_base(but_base,/row,/exclusive,/align_center)
slis_base=widget_base(but_base,/row,/align_center)
lat_lon_but=widget_button(buts_base,value='Lat & Lon',uvalue='LAT_LON',/no_release)
alt_lon_but=widget_button(buts_base,value='Alt & Lon',uvalue='ALT_LON',/no_release)
alt_lat_but=widget_button(buts_base,value='ALT & LAT',uvalue='ALT_LAT',/no_release)
sli=widget_slider(slis_base,/drag,maximum=nlons,minimum=0,uvalue='SLI',xsize=260)

txt_base=widget_base(but_base,/row,/align_center)
lab=widget_label(txt_base,value='Location')
txt=''
slider_txt=widget_text(txt_base,value=txt,uvalue='LOC',xsize=12)

grid_base=widget_base(but_base,/row,/align_left)
lab=widget_label(grid_base,value='Show Grid ?       ')
grid_buts=widget_base(grid_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
grid_yes_but=widget_button(grid_buts,value='Yes',uvalue='GRIDYES',/no_release)
grid_no_but=widget_button(grid_buts,value='No',uvalue='GRIDNO',/no_release)

ghost_base = widget_base(but_base,/row,/align_left)
lab        = widget_label(ghost_base,value='Show Ghost Cells? ')
gh_buts    = widget_base(ghost_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_but    = widget_button(gh_buts,value='Yes',uvalue='YES',/no_release)
no_but     = widget_button(gh_buts,value='No',uvalue='NO',/no_release)

polar_base   = widget_base(but_base,/row,/align_left)
lab          = widget_label(polar_base,value='Polar Plot ?')
polar_buts   = widget_base(polar_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_polarbut = widget_button(polar_buts,value='Yes', uvalue='POLAR', /no_release)
no_polarbut  = widget_button(polar_buts,value='No',  uvalue='CART',  /no_release)

log_base   = widget_base(but_base,/row,/align_left)
lab        = widget_label(log_base,value='Plot Log (Log10) ?')
log_buts   = widget_base(log_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_logbut = widget_button(log_buts,value='Yes',uvalue='LOGYES',/no_release)
no_logbut  = widget_button(log_buts,value='No',uvalue='LOGNO',/no_release)


vector_base=widget_base(but_base,/row,/align_left)
lab=widget_label(vector_base,value='Plot Vectors ?    ')
vector_yesno_buts=widget_base(vector_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
vyes_but=widget_button(vector_yesno_buts,value='Yes',uvalue='VYES',/no_release)
vno_but=widget_button(vector_yesno_buts,value='No',uvalue='VNO',/no_release)
vector_buts=widget_base(vector_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
vi_but=widget_button(vector_buts,value='Vi',uvalue='VI',/no_release)
vn_but=widget_button(vector_buts,value='Vn',uvalue='VN',/no_release)

vector_sli_base=widget_base(but_base,/column,/align_left,GRID_LAYOUT=1)
vector_cf_base = widget_base(vector_sli_base,/row)
cf_lab=widget_label(vector_cf_base,value ='Calculate Factor: ')
cf_sli=widget_slider(vector_cf_base,/drag,maximum=10,minimum=1,uvalue='CF',xsize=200)
vector_vs_base = widget_base(vector_sli_base,/row)
vs_lab = widget_label(vector_vs_base,value='Vectors Step:     ')
vs_sli=widget_slider(vector_vs_base,/drag,maximum=5,minimum=1,uvalue='VS',xsize=200)

cursor_base=widget_base(but_base,/row,/align_left)
cur_but=widget_button(cursor_base,value='Set Cursor',uvalue='SETCURSOR',/no_release)
lab=widget_label(cursor_base,value='X:')
txt=''
curx_txt=widget_text(cursor_base,value=txt,uvalue='CURX',xsize=12)
lab=widget_label(cursor_base,value='Y:')
txt=''
cury_txt=widget_text(cursor_base,value=txt,uvalue='CURY',xsize=12)

third_base=widget_base(lists_base,/column)

txt_base=widget_base(third_base,/row,/align_center)
lab=widget_label(txt_base,value='Mini')
txt=''
min_txt=widget_text(txt_base,value=txt,uvalue='MINI',xsize=12,/editable)
lab=widget_label(txt_base,value='Maxi')
max_txt=widget_text(txt_base,value=txt,uvalue='MAXI',xsize=12,/editable)
widget_control,sli,sensitive=0

base=widget_base(third_base,/column)
buts=widget_base(base,/row)
pc_base=widget_base(base,/row,/align_right)
lab=widget_label(buts,value='Write to a file?')
yn_buts=widget_base(buts,/row,/exclusive)
nowrite_but=widget_button(yn_buts,value='No',uvalue='NOWRITE',/no_release)
nowrite_cnt=1
widget_control,buts,/set_button
widget_control,nowrite_but,sensitive=0
yeswrite_but=widget_button(yn_buts,value='Yes',uvalue='YESWRITE',/no_release)
yeswrite_cnt=0
txt=''
writefile_txt=widget_text(buts,value=txt,uvalue='WRITEFILE',xsize=36,/editable)
okID=widget_button(pc_base,value='PLOT',uvalue='PLOT',xsize=100)
cancelID=widget_button(pc_base,value='CANCEL',uvalue='CAN',xsize=100)

no_cnt=1
;widget_control,gh_buts,/set_button
widget_control,no_but,sensitive=0

no_logcnt=1
;widget_control,log_buts,/set_button
widget_control,no_logbut,sensitive=0
plot_cnt=0
cursor_cnt=0

;widget_control,polar_buts,/set_button
polar = 0

gridyes_cnt=0
;widget_control,grid_buts,/set_button
widget_control,grid_no_but,sensitive=1

vyes_cnt=0
vno_cnt=1
;widget_control,vector_yesno_buts,/set_button
widget_control,vno_but,sensitive=1
widget_control,vector_buts,sensitive=0

cf_cnt=0
vs_cnt=0
;widget_control,cf_sli,/set_button
widget_control,cf_sli,sensitive=0
;widget_control,vs_sli,/set_button
widget_control,vs_sli,sensitive=0

if nowrite_cnt eq 1 then begin
	widget_control,writefile_txt,sensitive=0
endif

ptr = ptr_new({status:0,lists_base:lists_base,file_txt:file_txt,rb_txt:rb_txt,cb_txt:cb_txt,$
	rb:'',cb:'',vars_list:vars_list,filename:'',vars:ptr_new(vars),data:ptr_new(data),$
	nvars:0,nalts:0,nlats:0,nlons:0,lat_lon_cnt:0,alt_lon_cnt:0,alt_lat_cnt:0,$
	sli:sli,lat_lon_but:lat_lon_but,alt_lon_but:alt_lon_but,lat_lon_max:0,$
	alt_lon_max:0,alt_lat_max:0,alt_lat_but:alt_lat_but,sel:-1,nfiles:0,$
	slider_txt   : slider_txt,   $
        min_txt      : min_txt,      $
        max_txt      : max_txt,      $
        gh_buts      : gh_buts,      $
        no_but       : no_but,       $
	yes_but      : yes_but,      $
        no_cnt       : 1,            $
        yes_cnt      : 0,            $
        yes_logbut   : yes_logbut,   $
        bl_cnt       : 0,            $
	no_logbut    : no_logbut,    $
	yes_polarbut : yes_polarbut, $
	no_polarbut  : no_polarbut,  $
        polar        : polar,        $
        no_logcnt:1,yes_logcnt:0,log_buts:log_buts,$
	plot_cnt:0,cursor_x:0,cursor_y:0,curx_txt:curx_txt,cury_txt:cury_txt,$
	cursor_cnt:0,nowrite_cnt:1,yeswrite_cnt:0,writefile_txt:writefile_txt,$
	yeswrite_but:yeswrite_but,nowrite_but:nowrite_but,yn_buts:yn_buts,$
	grid_yes_but:grid_yes_but,grid_no_but:grid_no_but,vyes_but:vyes_but,vector_buts:vector_buts,$
	vno_but:vno_but,vi_but:vi_but,vn_but:vn_but,gridyes_cnt:0,vyes_cnt:0,$
	vi_cnt:0,vn_cnt:1,cf_sli:cf_sli,vs_sli:vs_sli,cf_cnt:0,vs_cnt:0,vs:1,cf:1})

widget_control,mainbase,set_uvalue=ptr
widget_control,mainbase,/realize
Xmanager,'thermosphere_vector2',mainbase
ptr_free,ptr
end

