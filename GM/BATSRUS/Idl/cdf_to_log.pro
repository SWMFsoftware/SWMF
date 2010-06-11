pro cdf_to_log,cdffilename,logfilename

; Convert a CDF file into an ASCII LOG file with date, position, B field

close,1
openw,1,logfilename
printf,1,"Extracted from ",cdffilename
printf,1,"year mo dy hr mn sc msc X Y Z bx by bz #START"

str2arr,cdffilename,cdffilenames,ncdf
for icdf=0,ncdf-1 do begin
    id = cdf_open(cdffilenames(icdf))
    re = cdf_inquire(id)
    nvars = re.nvars
    nrows = re.maxrec
    zvars = 0

    if (nrows lt 1) then begin
        cdf_control, id, variable=0, /zvariable, get_var_info=v
        zvars = 1
        nrows = v.maxrec
        if nrows eq 0 then nrows = v.maxrecs
        nvars = re.nzvars
    endif

    print, "number of vars and rows=",nvars, nrows

    if nvars lt 1 then retall
    if nrows lt 1 then retall

    timevar = "unknown"
    posvar  = "unknown"
    bvar    = "unknown"

    for i=0,nvars-1 do begin
        variables = cdf_varinq(id, i, zvar=zvars)
        print, i,'. ',variables.name
        
        if variables.name eq "IB_vector" or $
          variables.name eq "B_GSM_c"   or $
          variables.name eq "B_GSM"     or $
          variables.name eq "BGSMc"     or $
          variables.name eq "BGSM"      or $
          strpos(variables.name, "B_xyz_gse") eq 0 $
          then bvar = variables.name

        if timevar eq 'unknown' and strpos(variables.name,"Epoch") gt -1 $
          then timevar = variables.name

        if variables.name eq "SC_pos_GSM" or $
          variables.name eq "SC_pos_sm" or $
          variables.name eq "POS_GSM" or $
          variables.name eq "PGSM" or $
          variables.name eq "POS" $
          then posvar = variables.name

    endfor

;;    print,'reading b0 from polar file !!!'
;;    bvar = "MBCIGRF_GSM"

    if bvar eq "unknown" then begin
        print,"Could not find B variable in cdf file"
        retall
    endif

    if posvar eq "unknown" then begin
        print,"Could not find position variable in cdf file"
        retall
    endif

    if timevar eq "unknown" then begin
        print,"Could not find time variable in cdf file"
        retall
    endif

    cdf_varget, id, timevar, cdftime, rec_count=nrows
    cdf_varget, id, bvar, cdfbfield, rec_count=nrows
    cdf_varget, id, posvar, cdfpos, rec_count=nrows

    help,timevar,cdftime,posvar,cdfpos,bvar,cdfbfield

; Make sure timevar is a 1D array or make it one
    siz = size(cdftime)
    if siz(0) eq 2 then cdftime = reform(cdftime(0,*))
    
; Work around endian bug of IDL 7.0 if necessary
    if max(cdftime) lt 1.0 then begin
       cdftime   = swap_endian(cdftime)
       cdfbfield = swap_endian(cdfbfield)
       cdfpos    = swap_endian(cdfpos)
    endif

    if strpos(cdffilenames(icdf),'ge') gt -1 then begin
        print,"Geotail ? Dividing B by 10 !"
        cdfbfield = 0.1*cdfbfield
    endif

; Convert to Re units except for wind and polar satellites
    if strpos(cdffilenames(icdf),'wi_') lt 0 and $
      strpos(cdffilenames(icdf),'po_') lt 0 then begin
        print,"Converting position from km to Re"
        cdfpos = cdfpos/6378.0
    endif

; Convert to GSM if necessary
    if strpos(strlowcase(bvar)  , 'sm') lt 0 then begin
        print,'Converting bfield from GSE to GSM'
        gse_gsm, cdfbfield, cdftime
    endif
    if strpos(strlowcase(posvar), 'sm') lt 0 then begin
        print,'Converting position from GSE to GSM'
        gse_gsm, cdfpos, cdftime
    endif

    ; Print out data
    for i=0,nrows-1 do begin
;        if max(abs(cdfbfield(*,i))) lt 1000.0 then begin
            cdf_epoch, cdftime(i),year,month,day,hour,minute,second,milli,/break
            printf,1,$
              year,month,day,hour,minute,second,milli, $
              cdfpos(*,i),                             $
              cdfbfield(*,i),                          $
              format='(i5,5i3,i4,3f13.5,3e13.5)'
;        endif
    endfor
endfor

close,1
print,'Finished writing ',logfilename

end
