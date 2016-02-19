;For reading magnetograms in FITS or BATS-R-US format.

FUNCTION read_magnetogram, file, PlotRadius, UseBATS

  if UseBATS then begin

     gettype, file, filetype, npictinfile
     openfile, 10, file, filetype
     get_pict, 10, file, filetype, 1, x, var, error
     close, 10

     if gencoord then begin
        print, 'file '+file+' should contain a regular grid'
        retall
     endif

     case ndim of
        2:begin
           if PlotRadius ne 1.0 then begin
              print,'PlotRadius cannot be specified with 2D data!'
              retall
           endif

           if variables(0) ne "Longitude" or variables(1) ne "Latitude" or $
              variables(2) ne 'Br' then begin
              print, 'variables should be Longitude Latitude Br!'
              retall
           endif

           nlon = nx[0]
           nlat = nx[1]
           
           mag_info = {nlon:nlon,$
                       nlat:nlat,$
                       longitude:fltarr(nlon,nlat),$
                       latitude:fltarr(nlon,nlat),$
                       br_field:fltarr(nlon,nlat),$
                       bphi_field:fltarr(nlon,nlat),$
                       btheta_field:fltarr(nlon,nlat)}

           mag_info.longitude = x(*,*,0)*!dtor
           mag_info.latitude  = x(*,*,1)*!dtor
           mag_info.br_field = var(*,*,0)           
           
        end

        3:begin
           if variables(0) ne "Radius" or variables(1) ne "Longitude" or $
              variables(2) ne "Latitude" or variables(3) ne 'Br' then begin
              print, 'variables should be Radius Longitude Latitude Br!'
              retall
           endif
           
           nlon = nx[1] - 1
           nlat = nx[2]
           
           mag_info = {nlon:nlon,$
                       nlat:nlat,$
                       longitude:fltarr(nlon,nlat),$
                       latitude:fltarr(nlon,nlat),$
                       br_field:fltarr(nlon,nlat),$
                       bphi_field:fltarr(nlon,nlat),$
                       btheta_field:fltarr(nlon,nlat)}

           radius = x(*,0,0,0)
           longitude = x(0,*,*,1)
           latitude  = x(0,*,*,2)
           
           ; find index for the cut                                             
           d = abs(radius - PlotRadius)
           icut = where( d eq min(d) )
           br_field     = var(icut,*,*,0)
           bphi_field   = var(icut,*,*,1)
           btheta_field = var(icut,*,*,2)
          
           mag_info.br_field = reform(br_field[0,0:nlon-1,*])
           mag_info.bphi_field = reform(bphi_field[0,0:nlon-1,*])
           mag_info.btheta_field = reform(btheta_field[0,0:nlon-1,*])
           mag_info.longitude = reform(longitude[0,0:nlon-1,*])
           mag_info.latitude = reform(latitude[0,0:nlon-1,*])
        end
        else: begin
           print, 'ndim=', ndim, ' should be 2 or 3'
           retall
        end
     endcase
  endif else begin
                                ;read magnetogram in FITS format. For
                                ;transfering to Python, read_fits
                                ;function can be replaced by astropy
                                ;function.
     br_field=read_fits(file,index,/noscale)
     s=size(br_field)
     nlon=s[1]
     nlat=s[2]
     
     mag_info = {nlon:nlon,$
                       nlat:nlat,$
                       longitude:fltarr(nlon,nlat),$
                       latitude:fltarr(nlon,nlat),$
                       br_field:fltarr(nlon,nlat),$
                       bphi_field:fltarr(nlon,nlat),$
                       btheta_field:fltarr(nlon,nlat)}

     lat=findgen(nlat)*2./nlat
     lat=asin(lat-lat[nlat-1]/2.)
     lon=findgen(nlon)*!DPI*2./nlon
     latitude=fltarr(nlon,nlat)
     longitude=fltarr(nlon,nlat)
     for i=0,nlon-1 do begin
        for j=0,nlat-1 do begin
           latitude[i,j]=lat[j]
           longitude[i,j]=lon[i]
        endfor
     endfor
     
     mag_info.longitude=longitude
     mag_info.latitude=latitude
     mag_info.br_field=br_field
  endelse

return, mag_info

end
