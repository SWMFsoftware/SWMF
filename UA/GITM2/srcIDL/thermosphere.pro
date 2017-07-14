
pro read_thermosphere_file, filelist, nvars, nalts, nlats, nlons,vars,data

  f = filelist(0)
  all = findfile('*'+strmid(filelist(0),6,18))

  nfiles = n_elements(all)

  if (nfiles eq 16) then begin
    nBLKlat = 4
    nBLKlon = 4
  endif

  if (nfiles eq 2) then begin
    nBLKlat = 1
    nBLKlon = 2
  endif

  if (nfiles eq 32) then begin
    nBLKlat = 4
    nBLKlon = 8
  endif

  for n=0,nfiles-1 do begin

      file = all(n)

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

      help, data

      for k = 0, nalts-1 do begin
          for j = 0, nlats-1 do begin
              for i = 0, nlons-1 do begin
                  readf,1,tmp, format=format
                  ii = (n mod nBLKlon)*(nlons-4) + i
                  jj = (n / nBLKlon)*(nlats-4) + j
                  data(*,ii,jj,k) = tmp
              endfor
          endfor	
      endfor
      close,1
  endfor

  nlons = nlo
  nlats = nla

end

;Friday 27 September 2002
;Programmer: Ying Luo
;This program has widget interface and can do the three dimensional types plot:
;	1. Lat/Long
;	2. Alt/Long
;	3. Alt/Lat
;	4. Select a variable from 63 variables.
;	5. Show ghost cells
;	6. Plot Log?
;	7. Show Grid when display the images.
;	8. Set cursor's position and do plot according to the position.
;	9. Use position(POS=[0.05,0.05,0.8,0.96]) in the Contour procedure.
;	10.A color bar beside the contour picture.
;	11.Write the result to a file or not? If user selected to write to a 
;          file,after do plot again, program will write the result to the file.

pro thermosphere_event,event

widget_control,event.top,get_uvalue=ptr,/no_copy
widget_control,event.id,get_uvalue=whichevent

case whichevent of
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
                      read_thermosphere_file, filelist, nvars, nalts, nlats, nlons,vars,data
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
		  endelse
		endelse	
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
		datatoplot=reform(data(sel,*,*,k))
		maxi=max(datatoplot)
		mini=min(datatoplot)
		x=reform(lon(*,*,k))
		y=reform(lat(*,*,k))
                location = string(alt(0,0,k),format='(f5.1)')+' km Altitude'
		xtitle='Longitude (deg)'
		ytitle='Latitude (deg)'
		if yes eq 1 then begin
		   xrange=mm(lon)
		   yrange=mm(lat)
		endif
		if yes eq 0 then begin
		   xrange=[0,360]
		   yrange=[-90,90]
		endif
		ygrids=nlats
		xgrids=nlons
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
;	   readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
           makect,'mid'
	   clevels = findgen(31)/30 * 253.0 + 1.0

	   ;Draw color bar.
	   pos = [0.82,0.05,0.87,0.96]
	   maxmin = mm(levels)
	   title = ''
	   plotct,254,pos,maxmin,title,/right,color=color

	   ;If user DOES NOT do the Ghost Cells Contour. 
	   if yes eq 0 then begin
	   ;If user DO NOT do plot log. 
	   if nolog eq 1 then begin
		maxi=max(datatoplot)
		mini=min(datatoplot)
		if mini eq maxi then maxi=mini*1.01+1.0
	   	levels = findgen(31)/30.0*(maxi-mini) + mini

	      contour,datatoplot(*,*),x(*,*),y(*,*),POS=[0.08,0.08,0.8,0.96],$
		levels=levels,xstyle=1,ystyle=1,$
		xrange=xrange,yrange=yrange,$
		title=variable+' at '+location,c_colors=clevels,$
		xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE

	    	;Draw grid.
	    	for i=0,ygrids-1 do begin
 			oplot,mm(x),[y(0,i),y(0,i)]
	    	endfor
	    	for j=0,xgrids-1 do begin
			oplot,[x(j,0),x(j,0)],mm(y)
	    	endfor
		
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
	     endif   ;End if nolog eq 0.

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

		contour,datatoplot(*,*),x(*,*),y(*,*),POS=[0.08,0.08,0.8,0.96],$
		  levels=levels,xstyle=1,ystyle=1,$
		  xrange=xrange,yrange = yrange,$
		  title=variable+' at '+location,c_colors=clevels,$
		  xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE

		;Draw grid.
	    	for i=0,ygrids-1 do begin
 			oplot,mm(x),[y(0,i),y(0,i)]
	    	endfor
	    	for j=0,xgrids-1 do begin
			oplot,[x(j,0),x(j,0)],mm(y)
	    	endfor
		
	      	print,'Not Show Ghost Cells and do plot log'
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
	     endif	;End if yeslog eq 1
	   endif 	;ENd if yes eq 0

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

	        contour,datatoplot(*,*),x(*,*),y(*,*),POS=[0.08,0.08,0.8,0.96],$
		   levels=levels,xstyle=1,ystyle=1,$
		   xrange=xrange,yrange=yrange,$
		   title='Contour Plot Thermosphere',c_colors=clevels,$
		   xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE
	
		;Draw grid.
	    	for i=0,ygrids-1 do begin
 			oplot,mm(x),[y(0,i),y(0,i)]
	    	endfor
	    	for j=0,xgrids-1 do begin
			oplot,[x(j,0),x(j,0)],mm(y)
	    	endfor
		
		print,'Show Ghost Cells and no plot log'
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
	     endif	;End of if nolog eq 1

	     ;If user does want to do the Plot Log. 
	     if yeslog eq 1 then begin
		nzsubs=where(datatoplot gt 0, cnt)
		datatoplot(nzsubs)=datatoplot(nzsubs)
		datatoplot=ALOG10(datatoplot)
		maxi=max(datatoplot)
		mini=min(datatoplot)
		if mini eq maxi then maxi=mini+1
	   	levels = findgen(31)/30.0*(maxi-mini) + mini	
			
		contour,datatoplot(*,*),x(*,*),y(*,*),POS=[0.08,0.08,0.8,0.96],$
		levels=levels,xstyle=1,ystyle=1,$
		xrange=xrange,yrange=yrange,$
		title='Contour Plot Thermosphere',c_colors=clevels,$
		xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE

		;Draw grid.
	    	for i=0,ygrids-1 do begin
 			oplot,mm(x),[y(0,i),y(0,i)]
	    	endfor
	    	for j=0,xgrids-1 do begin
			oplot,[x(j,0),x(j,0)],mm(y)
	    	endfor
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
	     endif	;End of if yeslog eq 1
	   endif	;End of if yes eq 1

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

pro thermosphere

title='Plot Thermoshpere'
mainbase=widget_base(title=title,/column)
base = widget_base(mainbase,/column)

dirlab=widget_base(base,/column,/frame)
lab=widget_label(dirlab,value='Enter file name (can use *) [p0000b0001i0200.dat]:',/align_left)
filename_base=widget_base(dirlab,/row)
filelist = findfile("-t *.dat")
txt=filelist(0)
file_txt=widget_text(filename_base,value=txt,uvalue='FILENAME',xsize=46,/editable,/align_left)
sear_but=widget_button(filename_base,value='Search',uvalue='SEARCH',xsize=100)

lists_base=widget_base(base,/row,/frame)
vars_base=widget_base(lists_base,/column)
sel=-1
lab=widget_label(vars_base,value='Variables',/align_center)
vars=''
vars_list=widget_list(vars_base,value=vars,uvalue='SELVARS',xsize=20,$
	ysize=80,/align_right)
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


ghost_base=widget_base(but_base,/row,/align_left)
lab=widget_label(ghost_base,value='Show Ghost Cells?')
gh_buts=widget_base(ghost_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_but=widget_button(gh_buts,value='Yes',uvalue='YES',/no_release)
no_but=widget_button(gh_buts,value='No',uvalue='NO',/no_release)

log_base=widget_base(but_base,/row,/align_left)
lab=widget_label(log_base,value='Plot Log (Log10) ?')
log_buts=widget_base(log_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_logbut=widget_button(log_buts,value='Yes',uvalue='LOGYES',/no_release)
no_logbut=widget_button(log_buts,value='No',uvalue='LOGNO',/no_release)

cursor_base=widget_base(but_base,/row,/align_left)
cur_but=widget_button(cursor_base,value='Set Cursor',uvalue='SETCURSOR',/no_release)
lab=widget_label(cursor_base,value='X:')
txt=''
curx_txt=widget_text(cursor_base,value=txt,uvalue='CURX',xsize=12)
lab=widget_label(cursor_base,value='Y:')
txt=''
cury_txt=widget_text(cursor_base,value=txt,uvalue='CURY',xsize=12)

txt_base=widget_base(butt_base,/row,/align_center)
lab=widget_label(txt_base,value='Mini')
txt=''
min_txt=widget_text(txt_base,value=txt,uvalue='MINI',xsize=12)
lab=widget_label(txt_base,value='Maxi')
max_txt=widget_text(txt_base,value=txt,uvalue='MAXI',xsize=12)
widget_control,sli,sensitive=0

base=widget_base(mainbase,/column)
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
widget_control,gh_buts,/set_button
widget_control,no_but,sensitive=0

no_logcnt=1
widget_control,log_buts,/set_button
widget_control,no_logbut,sensitive=0
plot_cnt=0
cursor_cnt=0

if nowrite_cnt eq 1 then begin
	widget_control,writefile_txt,sensitive=0
endif

ptr = ptr_new({status:0,lists_base:lists_base,file_txt:file_txt,$
	vars_list:vars_list,filename:'',vars:ptr_new(vars),data:ptr_new(data),$
	nvars:0,nalts:0,nlats:0,nlons:0,lat_lon_cnt:0,alt_lon_cnt:0,alt_lat_cnt:0,$
	sli:sli,lat_lon_but:lat_lon_but,alt_lon_but:alt_lon_but,lat_lon_max:0,$
	alt_lon_max:0,alt_lat_max:0,alt_lat_but:alt_lat_but,sel:-1,nfiles:0,$
	slider_txt:slider_txt,min_txt:min_txt,max_txt:max_txt,gh_buts:gh_buts,no_but:no_but,$
	yes_but:yes_but,no_cnt:1,yes_cnt:0,yes_logbut:yes_logbut,$
	no_logbut:no_logbut,no_logcnt:1,yes_logcnt:0,log_buts:log_buts,$
	plot_cnt:0,cursor_x:0,cursor_y:0,curx_txt:curx_txt,cury_txt:cury_txt,$
	cursor_cnt:0,nowrite_cnt:1,yeswrite_cnt:0,writefile_txt:writefile_txt,$
	yeswrite_but:yeswrite_but,nowrite_but:nowrite_but,yn_buts:yn_buts})

widget_control,mainbase,set_uvalue=ptr
widget_control,mainbase,/realize
Xmanager,'thermosphere',mainbase
ptr_free,ptr
end






























