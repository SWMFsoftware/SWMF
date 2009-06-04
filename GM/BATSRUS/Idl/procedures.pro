;^CFG COPYRIGHT VAC_UM
; Written by G. Toth for the Versatile Advection Code and BATSRUS
; Some improvements by Aaron Ridley.
;
; Procedures for
;
; reading ascii and binary data produced by VAC, VACINI, BATSRUS etc:
;    openfile, gettype, gethead, get_pict, 
;    get_pict_asc, get_pict_bin, get_pict_log, get_log
; reading numbers and strings from input:
;    asknum, askstr, str2arr, arr2arr, readplotpar, readlimits
; transforming initial data:
;    regulargrid, polargrid, unpolargrid, spheregrid, getaxes
;    interpol_logfiles, interpol_log
; calculating functions of the data
;    getfunc, getlimits
; plotting
;    plot_func, plotgrid, plot_log
; calculating cell corners and cell volumes for general 2D grids
;    gengrid
; comparing two w,x or wlog arrays for relative differences
;    compare, rms_logfiles, rel_error
; checking index ranges for functions quadruplet and triplet
;    checkdim
; procedure "quit" as an alias for "exit"
;    quit
;
; Functions for
;
; calculating first derivative in 1D
;    diff1
; calculating derivatives in 2D for Cartesian grids to 2nd,3rd,4th order
;    diff2,diff3,diff4
; calculating the laplace on 1D or 2D Cartesian grids
;    laplace
; calculating minmod limited slope
;    minmod
; calculating symmetric differences with respect to some mirror plane
;    symmdiff
; calculating derivatives in 2D for general grids
;    grad,div,curl,grad_rz,div_rz,curl_rz, filledge,intedge,intedge_rz
; taking a part of an array or coarsen an array
;    triplet, quadruplet, coarsen
; eliminating degenerate dimensions from an array
;    reform2
; converting logfile time or date+time into hours, getting other functions
;    log_time, log_func


;=============================================================================
function reform2,x

;Remove all degenerate dimensions from x

if n_elements(x) lt 2 then return,x

siz=size(x)
siz=siz(1:siz(0))
return,reform(x,siz(where(siz gt 1)))

end

;=============================================================================
function log_time,wlog,wlognames,timeunit

; Obtain time in hours from wlog and wlognames
; If the log file contains 't' or 'time', simply convert seconds to hours
; If the log file contains date and time, calculate time relative to
; the beginning of the first day. 
; This algorithm works only if the whole file is in the same month.

nwlog = n_elements(wlognames)
hours = 0.0*wlog(*,0)

if nwlog lt 0 then return, hours

istep = -1
iyear = -1
itime = -1
ihour = -1
for i = 0, nwlog-1 do begin
    if wlognames(i) eq 'step' or wlognames(i) eq 'it' then istep=i
    if wlognames(i) eq 'year' or wlognames(i) eq 'yr' then iyear=i
    if wlognames(i) eq 'time' or wlognames(i) eq 't'  then itime=i
    if wlognames(i) eq 'hour' or wlognames(i) eq 'hours'  then ihour=i
endfor

if itime gt -1 then begin
    hours = wlog(*,itime)/3600.0
endif else begin
    if iyear eq -1 then begin
        if ihour gt -1 then begin
            hours = wlog(*,ihour)
            for i=1, n_elements(hours)-1 do $
              while hours(i) lt hours(i-1) do $
              hours(i) = hours(i) + 24.0
        endif
    endif else begin
        iday  = iyear+2
        ihour = iday+1
        imin  = iday+2
        isec  = iday+3
        imsc  = iday+4
        hours = (wlog(*,iday)-wlog(0,iday))*24.0 + $
          wlog(*,ihour) + wlog(*,imin)/60.0 + wlog(*,isec)/3600.0 + $
          wlog(*,imsc)/3.6e6
    endelse
endelse

if max(hours) eq min(hours) then begin
    if istep gt -1 then begin
        print, 'Could not find time information, using steps'
        hours = wlog(*,istep)
    endif else begin
        print, 'Could not find time information, using line number'
        hours = findgen(n_elements(wlog(*,0)))
    endelse
endif

logtime = hours

if n_elements(timeunit) gt 0 then begin
    case timeunit of
        '1': logtime = hours*3600
        's': logtime = hours*3600
        'm': logtime = hours*60
        else: 
    endcase
endif

return, logtime
end

;=============================================================================
function log_func, wlog, varnames, varname, error

    error = 0
    ivar  = where(varnames eq varname) & ivar = ivar(0)
    ; Variable is found, return with array
    if ivar ge 0 then return, wlog(*,ivar)

    ; Try calculating temperature or pressure
    if varname eq 'T' then begin
                                ; Convert p[nPa]/n[/cc] to T[eV]
        ivar = where( varnames eq 'p')   & ivar = ivar(0)
        jvar = where( varnames eq 'rho') & jvar = jvar(0)
        if ivar ge 0 and jvar ge 0 then $
          return,6241.5*wlog(*,ivar)/wlog(*,jvar)

    endif else if varname eq 'p' then begin
                                ; Convert T[eV]*n[/cc] to p[nPa] to T[eV]
        ivar = where( varnames eq 'T')   & ivar = ivar(0)
        jvar = where( varnames eq 'rho') & jvar = jvar(0)
        if ivar ge 0 and jvar ge 0 then $
          return, wlog(*,ivar)*wlog(*,jvar)/6241.5
    endif

    error = 1
    return,0*wlog(*,0)
    
end

;=============================================================================
pro openfile,unit,filename,filetype

   on_error,2

   close,unit
   case filetype of
       'log'   :openr,unit,filename
       'ascii' :openr,unit,filename
       'binary':openr,unit,filename,/f77_unf
       'real4' :openr,unit,filename,/f77_unf
       'BINARY':openr,unit,filename,/f77_unf
       'REAL4' :openr,unit,filename,/f77_unf
       else    :print,'Openfile: unknown filetype:',filetype
   endcase
end

;=============================================================================
pro gettype,filenames,filetypes,npictinfiles

on_error,2

filetypes=filenames
npictinfiles=intarr(n_elements(filenames))
for ifile=0,n_elements(filenames)-1 do begin
    if   strpos(filenames(ifile),'.log') gt 0 $
      or strpos(filenames(ifile),'.sat') gt 0 then begin
        filetypes(ifile)    = 'log'
        npictinfiles(ifile) = 1
    endif else begin
        ; Obtain filetype based on the length info in the first 4 bytes
        close,10
        openr,10,filenames(ifile)
        lenhead=long(1)
        readu,10,lenhead
        if lenhead ne 79 and lenhead ne 500 then ftype='ascii' else begin
            ; The length of the 2nd line decides between real4 and binary
            ; since it contains the time, which is real*8 or real*4
            head=bytarr(lenhead+4)
            len=long(1)
            readu,10,head,len
            case len of
                20: ftype='real4'
                24: ftype='binary'
                else: begin
                    print,'Error in GetType: strange unformatted file:',$
                      filenames(ifile)
                    retall
                end
            endcase
            if lenhead eq 500 then ftype = strupcase(ftype)
        endelse
        close,10
        
        ; Obtain file size and number of snapshots
        openfile,1,filenames(ifile),ftype
        status=fstat(1)
        fsize=status.size

        pointer=0
        pictsize=1
        npict=0
        while pointer lt fsize do begin
            ; Obtain size of a single snapshot
            point_lun,1,pointer
            gethead,1,ftype,pictsize=pictsize
            npict=npict+1
            pointer=long64(pointer) + pictsize
        endwhile
        close,1

        npictinfiles(ifile)=npict
        filetypes(ifile)   =ftype
    endelse
  endfor
end

;=============================================================================
pro show_head,filename

on_error,2
if n_elements(filename) eq 0 then filename='file.out'
gettype,filename,filetype
openfile,1,filename,filetype
gethead,1,filetype,$
  headline,phys,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
close,1
print,'filename  = ',filename, format="(a,a)"
print,'filetype  = ',filetype, format="(a,a)"
print,'headline  = ',headline, format="(a,a)"
print,'phys      = ',phys,     format="(a,a)"
print,'it        = ',it,       format="(a,i8)"
print,'time      = ',time,     format="(a,g15.8)"
print,'gencoord  = ',gencoord, format="(a,i8)"
print,'ndim      = ',ndim,     format="(a,i8)"
print,'neqpar    = ',neqpar,   format="(a,i8)"
print,'nw        = ',nw,       format="(a,i8)"
print,'nx        = ',nx,       format="(a,3i8)"
print,'eqpar     = ',eqpar,    format="(a,100g15.8)"
print,'coords    =',variables(0:ndim-1)
print,'plotvars  =',variables(ndim:ndim+nw-1)
print,'eqpars    =',variables(ndim+nw:*)

end
;=============================================================================
pro gethead,unit,filetype,headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nx,eqpar,variables,pictsize=pictsize

;; on_error,2

ftype = strlowcase(filetype)

if ftype eq filetype then lenstr = 79 else lenstr = 500

;Type definitions
headline=''
for i=1, lenstr do headline=headline+' '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname=''
for i=1, lenstr do varname=varname+' '

;Remember pointer position at beginning of header
point_lun,-unit,pointer0

;Read header
case ftype of
    'log': begin
        readf,unit,headline
        headline=headline+'_mhd13'
        readf,unit,varname
        nw=n_elements(strsplit(varname))-2
        varname ='hour '+varname
        ; reset pointer
        point_lun,unit,pointer0
        it=0
        time=0.0
        gencoord=0
        neqpar=0
        eqpar=0.0
        ndim=1
        nx=lonarr(1)
        nx(0)=1
    end
    'ascii': begin
        time=double(1)
        readf,unit,headline
        readf,unit,it,time,ndim,neqpar,nw
        gencoord=(ndim lt 0)
        ndim=abs(ndim)
        nx=lonarr(ndim)
        readf,unit,nx
        eqpar=dblarr(neqpar)
        readf,unit,eqpar
        readf,unit,varname
    end
    'binary':begin
        time=double(1)
        readu,unit,headline
        readu,unit,it,time,ndim,neqpar,nw
        gencoord=(ndim lt 0)
        ndim=abs(ndim)
        nx=lonarr(ndim)
        readu,unit,nx
        eqpar=dblarr(neqpar)
        readu,unit,eqpar
        readu,unit,varname
    end
    'real4': begin
        time=float(1)
        readu,unit,headline
        readu,unit,it,time,ndim,neqpar,nw
        gencoord=(ndim lt 0)
        ndim=abs(ndim)
        nx=lonarr(ndim)
        readu,unit,nx
        eqpar=fltarr(neqpar)
        readu,unit,eqpar
        readu,unit,varname
    end
    else: begin
        print,'Gethead: unknown filetype',filetype
        retall
    end
endcase

if keyword_set(pictsize) then begin
                                ; Calculate the picture size
                                ; Header length
    point_lun,-unit,pointer1
    headlen=pointer1-pointer0
                                ; Number of cells
    nxs=1
    for idim=1,ndim do nxs=nxs*nx(idim-1)
                                ; Snapshot size = header + data + recordmarks
    case ftype of
        'log'   :pictsize = 1
        'ascii' :pictsize = headlen + (18*(ndim+nw)+1)*nxs
        'binary':pictsize = headlen + 8*(1+nw)+8*(ndim+nw)*nxs
        'real4' :pictsize = headlen + 8*(1+nw)+4*(ndim+nw)*nxs
    endcase
endif else begin
                                ; Get variables and physics
    variables=strsplit(varname,/extract)
    tmp=strsplit(strtrim(headline,2),'_',/extract)
    if n_elements(tmp) eq 2 then begin
        headline=tmp(0)
        physics=tmp(1)
                                ; work around for a bug
        if physics eq '123' or physics eq '223' then physics='mhd23'
    endif
endelse

end

;=============================================================================
pro get_pict,unit,filetype,npict,x,w,$
    headline,physics,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
    rBody,error

   on_error,2

   error=0

   if(eof(unit))then begin
      error=1
      return
   endif

   ; Get current pointer position
   point_lun,-unit,pointer

   ; Skip npict-1 snapshots
   ipict=0
   pictsize=1
   while ipict lt npict-1 and not eof(unit) do begin
      ipict=ipict+1
      gethead,unit,filetype,pictsize=pictsize
      pointer=long64(pointer) + pictsize
      point_lun,unit,pointer
   endwhile

   ; Backup 1 snapshot if end of file
   if eof(unit) then begin
       error=1
       point_lun,unit,pointer-pictsize
   endif

   ; Read header information
   gethead,unit,filetype,headline,physics,$
       it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables

   ; set rBody if listed among the parameters
   for i = nDim + nW, n_elements(variables)-1 do begin
       iPar = i - nDim - nW
       if variables(i) eq 'rbody' or variables(i) eq 'rBody' then $
           rBody = eqpar(iPar)
   endfor

   ; Read data
   case strlowcase(filetype) of
       'log':    get_pict_log ,unit, npict, ndim, nw, nx, x, w
       'ascii':  get_pict_asc ,unit, npict, ndim, nw, nx, x, w
       'binary': get_pict_bin ,unit, npict, ndim, nw, nx, x, w
       'real4':  get_pict_real,unit, npict, ndim, nw, nx, x, w
       else:    begin
           print,'get_pict: unknown filetype:',filetype
           error=1
           close,unit
       end
   endcase

   if ndim eq 2 and nx(ndim-1) eq 1 then begin
       ; sort x and w according to the x + factor * y function
       ; where factor is a transcendent number to avoid coincidences
       factor = exp(1.0d0)
       isort = sort( x(*,*,0) + factor*x(*,*,1) )
       x(*,0,*) = x(isort,0,*)
       w(*,0,*) = w(isort,0,*)
   endif

   if ndim eq 1 and gencoord then begin
       ; sort points to have increasing coordinates
       isort = sort(x)
       x     = x(isort)
       w     = w(isort,*)
   endif

end

;=============================================================================
pro get_pict_log,unit,npict,ndim,nw,nx,x,w

get_log, unit, w, wlognames, x, 's'
ndim = 1
nx(0)= n_elements(x)
nw   = n_elements(wlognames)

end

;=============================================================================
pro get_pict_asc,unit,npict,ndim,nw,nx,x,w

  on_error,2
  ;----------------------------------------
  ; Read coordinates and values row by row
  ;----------------------------------------
  wrow=dblarr(nw)
  xrow=dblarr(ndim)
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    x=dblarr(nx(0),ndim)
    w=dblarr(nx(0),nw)
    for x0=0L,nx(0)-1 do begin
      readf,unit,xrow,wrow
      x(x0,0:ndim-1)=xrow(0:ndim-1)
      w(x0,0:nw-1)  =wrow(0:nw-1)
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    x=dblarr(nx(0),nx(1),ndim)
    w=dblarr(nx(0),nx(1),nw)
    for x1=0L,nx(1)-1 do begin
      for x0=0L,nx(0)-1 do begin
        readf,unit,xrow,wrow
        x(x0,x1,0:ndim-1)=xrow(0:ndim-1)
        w(x0,x1,0:nw-1)  =wrow(0:nw-1)
      endfor
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    x=dblarr(nx(0),nx(1),nx(2),ndim)
    w=dblarr(nx(0),nx(1),nx(2),nw)
    for x2=0L,nx(2)-1 do begin
      for x1=0L,nx(1)-1 do begin
        for x0=0L,nx(0)-1 do begin
          readf,unit,xrow,wrow
          x(x0,x1,x2,0:ndim-1)=xrow(0:ndim-1)
          w(x0,x1,x2,0:nw-1)=wrow(0:nw-1)
        endfor
      endfor
    endfor
  end
  endcase
end

;=============================================================================
pro get_pict_bin,unit,npict,ndim,nw,nx,x,w

  on_error,2
  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=dblarr(n1,ndim)
    w=dblarr(n1,nw)
    wi=dblarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=dblarr(n1,n2,ndim)
    w=dblarr(n1,n2,nw)
    wi=dblarr(n1,n2)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=dblarr(n1,n2,n3,ndim)
    w=dblarr(n1,n2,n3,nw)
    wi=dblarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;=============================================================================
pro get_pict_real,unit,npict,ndim,nw,nx,x,w

  on_error,2
  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=fltarr(n1,ndim)
    w=fltarr(n1,nw)
    wi=fltarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=fltarr(n1,n2,ndim)
    w=fltarr(n1,n2,nw)
    readu,unit,x
    wi=fltarr(n1,n2)
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=fltarr(n1,n2,n3,ndim)
    w=fltarr(n1,n2,n3,nw)
    wi=fltarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;=============================================================================
pro asknum,prompt,var,doask

   on_error,2

   if var eq 0 then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+strtrim(string(var),2)+' ? ',tmp
         if tmp ne '' then reads,tmp,var
      endif else print,prompt,'=',var
   endelse
end

;=============================================================================
pro askstr,prompt,var,doask

   on_error,2

   if var eq '' then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+var+' ? ',tmp
         if tmp ne '' then var=tmp
      endif else print,prompt,'=',var
   endelse
end

;=============================================================================
pro str2arr,s,a,n,sep

; If s is an array copy it to a.
; If s is a string then split string s at the sep characters 
; (default is space) into array a.
; If n is 0, it will be the size of the array on output
; If n is defined, fill up the rest of the array with the last defined element
; If n is defined but smaller than the number of elements in s, print a warning

on_error,2

if n_elements(s) gt 1 then     $
    a0 = s                     $
else if keyword_set(sep) then  $
   a0=strsplit(s,sep,/extract) $
else                           $
   a0=strsplit(s,/extract)

n0=n_elements(a0)

if not keyword_set(n) then begin
   a=a0
   n=n0
endif else if n ge n0 then begin
   a=strarr(n)
   a(0:n0-1)=a0
   if n0 lt n then a(n0:n-1)=a0(n0-1)
endif else begin
   a=strarr(n)
   a=a0(0:n-1)
   print,'Warning: more than',n,' values defined by string: ',s,$
       FORMAT='(a,i3,a,a)'
endelse

end
;=============================================================================
pro arr2arr,a,n

; If a is a scalar or has fewer elements than n then 
; extend it to an array of n elements.
; If a has more elements than n, reduce the number of elements to n.

on_error,2

k = n_elements(a)
if      k gt n then a = a(0:n-1) $                 ; truncate array
else if k lt n then a = [a, a(k-1)+intarr(n-k) ]   ; extend array

end
;===========================================================================
pro readplotpar,ndim,cut,cut0,plotdim,nfunc,func,funcs,funcs1,funcs2,$
   nplot,plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask

   on_error,2

   ; Determine dimension of plots based on cut or ndim,
   ; calculate reduced cut0 array by eliminating degenerate dimensions
   if keyword_set(cut) then begin
      cut0=reform2(cut)
      siz=size(cut0)
      plotdim=siz(0)
   endif else begin
      plotdim=ndim
      cut0=0
   endelse

   askstr,'func(s) (e.g. rho p ux;uz bx+by -T) ',func,doask
   if plotdim eq 1 then begin
      print,'1D plotmode: plot'
      plotmode='plot'
   endif else begin
      if plotmode eq 'plot' then plotmode=''
      print,'2D plotmode: shade/surface/cont/tv/polar/velovect/vector/stream'
      print,'2D +options: bar,body,fill,grid,irr,label,log,mesh,over,white'
      askstr,'plotmode(s)                ',plotmode,doask
   endelse
   askstr,'plottitle(s) (e.g. B [G];J)',plottitle,doask
   askstr,'autorange(s) (y/n)         ',autorange,doask

   nfunc=0
   str2arr,func,funcs,nfunc
   str2arr,plotmode,plotmodes,nfunc
   str2arr,plottitle,plottitles,nfunc,';'
   str2arr,autorange,autoranges,nfunc

   nplot = nfunc
   funcs1=strarr(nfunc)
   funcs2=strarr(nfunc)
   for ifunc=0,nfunc-1 do begin
      func12=strsplit(funcs(ifunc),';',/extract)
      funcs1(ifunc)=func12(0)
      if n_elements(func12) eq 2 then funcs2(ifunc)=func12(1)
      if strpos(plotmodes(ifunc),'over') ge 0 then nplot=nplot-1
   endfor

end
;===========================================================================
pro readtransform,ndim,nx,gencoord,transform,nxreg,xreglimits,wregpad,$
    physics,nvector,vectors,grid,doask

   on_error,2

   if (gencoord or transform eq 'unpolar') and ndim eq 2 then begin
      if transform eq '' then begin
        transform='none'
        askstr,"transform (r=regular/p=polar/u=unpolar/n=none)",transform,1
      endif else $
        askstr,"transform (r=regular/p=polar/u=unpolar/n=none)",transform,doask

      ; Complete name
      case transform of
          'r': transform='regular'
          'p': transform='polar'
          'u': transform='unpolar'
          'n': transform='none'
         else:
      endcase
      ; Get transformation parameters and calculate grid
      case 1 of
        transform eq 'regular':begin
           print,'Generalized coordinates, dimensions for regular grid'
           if n_elements(nxreg) ne 2 then nxreg=[0,0]
           if n_elements(xreglimits) ne 4 then xreglimits=dblarr(4) $
           else xreglimits=double(xreglimits)
           nxreg0=nxreg(0)
           nxreg1=nxreg(1)
           asknum,'nxreg(0) (use negative sign to limit x)',nxreg0,doask
           if nxreg0 lt 0 then begin
               nxreg0=abs(nxreg0)
               xmin=0 & xmax=0
               asknum,'xreglimits(0) (xmin)',xmin,doask
               asknum,'xreglimits(2) (xmax)',xmax,doask
               xreglimits(0)=xmin
               xreglimits(2)=xmax
           endif
           asknum,'nxreg(1) (use negative sign to limit y)',nxreg1,doask
           if nxreg1 lt 0 then begin
               nxreg1=abs(nxreg1)
               ymin=0 & ymax=0
               asknum,'xreglimits(1) (ymin)',ymin,doask
               asknum,'xreglimits(3) (ymax)',ymax,doask
               xreglimits(1)=ymin
               xreglimits(3)=ymax
           endif
           grid=lindgen(nxreg0,nxreg1)
           nxreg=[nxreg0,nxreg1]
           wregpad=0
        end
        transform eq 'polar' or transform eq 'unpolar':begin
           getvectors,physics,nvector,vectors
           grid=lindgen(nx(0),nx(1))
        end
        transform eq 'none':grid=lindgen(nx(0),nx(1))
        else: print,'Unknown value for transform:',transform
      endcase
   endif else if (gencoord or transform eq 'unpolar') and ndim eq 3 then begin
       if transform eq '' then begin
           transform="none" & askstr, $
             "transform (p=polar/s=sphere/u=unpolar/n=none)",transform,1
       endif else $
         askstr, $
         "transform (p=polar/s=sphere/u=unpolar/n=none)",transform,doask
       case transform of
           'p': transform='polar'
           's': transform='sphere'
           'u': transform='unpolar'
           'n': transform='none'
           else:
       endcase
       case 1 of
           transform eq 'polar' or transform eq 'unpolar' or $
             transform eq 'sphere' :begin
               getvectors,physics,nvector,vectors
               grid=lindgen(nx(0),nx(1),nx(2))
           end
           transform eq 'none': grid=lindgen(nx(0),nx(1),nx(2))
           else: print,'Unknown value for transform:',transform
       endcase
   endif else case ndim of
       1: grid=lindgen(nx(0))
       2: grid=lindgen(nx(0),nx(1))
       3: grid=lindgen(nx(0),nx(1),nx(2))
   endcase

   ;===== GRID HELPS TO CREATE A CUT, E.G.: cut=grid(*,4)

   help,grid
end

;===========================================================================
pro getvectors,physics,nvector,vectors

   nvector = 0
   physic=strtrim(physics)
   if strlen(physics) lt 3 then return
   phys=strmid(physic,0,strlen(physic)-2)
   ndir=0
   reads,strmid(physic,strlen(physic)-1,1),ndir
   case phys of
   'rho':nvector=0
   'flx':nvector=0
   'hd' :begin
         nvector=1
         vectors=1
         end
   'hdadiab':begin
         nvector=1
         vectors=1
         end
   'mhdiso':begin
         nvector=2
         vectors=[1,ndir+1]
         end
   'mhd':begin
         nvector=2
         vectors=[1,ndir+2]
         end
   else:begin
      if nvector eq 0 then begin
         print,'Unrecognised physics: ',physics
         print,'Vector variables to transform for WREG'
         asknum,'nvector',nvector,doask
         if nvector gt 0 then begin
            vectors=intarr(nvector)
            read,PROMPT='Indices of first components in w? ',vectors
         endif
      endif
      end
   endcase
end

;===========================================================================
pro readlimits,nfunc,funcs,autoranges,noautorange,fmax,fmin,doask

   on_error,2

   if n_elements(fmax) ne nfunc then fmax=dblarr(nfunc) else fmax=double(fmax)
   if n_elements(fmin) ne nfunc then fmin=dblarr(nfunc) else fmin=double(fmin)

   ; check if there is any function for which autorange is 'y'
   noautorange=1
   for ifunc=0,nfunc-1 do noautorange=noautorange and autoranges(ifunc) eq 'n'

   if(noautorange)then begin
      for ifunc=0,nfunc-1 do begin
         f_min=fmin(ifunc)
         f_max=fmax(ifunc)
         asknum,'Min value for '+funcs(ifunc),f_min,doask
         asknum,'Max value for '+funcs(ifunc),f_max,doask
         fmin(ifunc)=f_min
         fmax(ifunc)=f_max
      endfor
   endif

end
;===========================================================================
pro getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,fmax,fmin,doask,$
                x,w,xreg,wreg,usereg,physics,eqpar,variables,cut,rcut

   on_error,2

   for ifunc=0,nfunc-1 do begin
      if autoranges(ifunc) eq 'n' then begin
         if first then begin
            f_min=fmin(ifunc)
            f_max=fmax(ifunc)
            asknum,'Min value for '+funcs(ifunc),f_min,doask
            asknum,'Max value for '+funcs(ifunc),f_max,doask
            fmin(ifunc)=f_min
            fmax(ifunc)=f_max
         endif
      endif else begin
         if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            xreg,wreg,physics,eqpar,variables,cut,rcut $
         else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            x,   w,   physics,eqpar,variables,cut,rcut

         f_max=max(f)
         f_min=min(f)
         if first then begin
            fmax(ifunc)=f_max
            fmin(ifunc)=f_min
         endif else begin
            if f_max gt fmax(ifunc) then fmax(ifunc)=f_max
            if f_min lt fmin(ifunc) then fmin(ifunc)=f_min
         endelse
      endelse
   endfor
end

;==============================================================================
function grid_data,x,y,data,nxreg,xreglimits,triangles,wregpad

return,griddata(x,y,data,$
        dimension=nxreg,$
	start=xreglimits(0:1),$
	delta=[(xreglimits(2)-xreglimits(0))/(nxreg(0)-1),  $
               (xreglimits(3)-xreglimits(1))/(nxreg(1)-1)] ,$
        triangles=triangles,$
;        method='NearestNeighbor',$
;        method='Linear',$
        method='InverseDistance',$
       smoothing=0.5,$
       max_per_sector=4, $
        missing=wregpad $
        )
end


;===========================================================================
pro regulargrid,x_old,nxreg_old,xreglimits_old,x,xreg,nxreg,xreglimits,$
                w,wreg,nw,wregpad,triangles,symmtri
;
;    Regularize grid and interpolate "w" via triangulate and trigrid.
;    The original "w" data is interpolated into "wreg", for points outside
;    the convex hull of the irregular grid the "wregpad(nw)" array is used.
;
;    If "x_old" and "x" or "nxreg_old" and "nxreg" are different
;    a triangulization is done first and a regular coordinate array
;    "xreg" is created. The size of the "xreg" array is in "nxreg(2)",
;    "xreglimits(4)" gives the limits. The triangles are saved in "triangles".
;
;    "q" can be interpolated from the irregular grid to the regular one by:
;
;    qreg(*,*)=trigrid(x(*,*,0),x(*,*,1),q,triangles,[0,0],xreglimits)

   on_error,2

   ;Floating underflow is not a real error, the message is suppressed
   err=check_math(1,1)

   xx=x(*,*,0)
   yy=x(*,*,1)

   ; Test distribution. If you uncomment the next lines you can
   ; take a look at the different "shock wave" representation
   ; on your grid for the 0-th variable (usually rho)
   ; for i=0L,n_elements(xx)-1 do $
   ;   if abs(xx(i))-0.2*abs(yy(i)) gt 0.01 then $
   ;       w(i,*,0)=2. else w(i,*,0)=1.

   ; Check if nxreg==nxreg_old and xreglimits==xreglimits_old and x==x_old
   newx=1
   nrectan=0
   if symmtri ne 1 and symmtri ne 2 then $
   if n_elements(nxreg_old) eq n_elements(nxreg) then $
   if max(abs(nxreg_old-nxreg)) eq 0 then $
   if n_elements(xreglimits) eq n_elements(xreglimits_old) then $
   if max(abs(xreglimits-xreglimits_old)) eq 0 then $
   if n_elements(x_old) eq n_elements(x) then $
   if max(abs(x_old-x)) eq 0 then newx=0

   if xreglimits(0) eq xreglimits(2) then begin
      xreglimits(0)=min(xx) & xreglimits(2)=max(xx)
   endif
   if xreglimits(1) eq xreglimits(3) then begin
      xreglimits(1)=min(yy) & xreglimits(3)=max(yy)
   endif

   if newx then begin
      print,'Triangulating...'
      x_old=x
      nxreg_old=nxreg
      xreglimits_old=xreglimits

      triangulate,float(xx),float(yy),triangles

      ; calculate conjugate triangulation and rectangles if required
      if symmtri eq 1 or symmtri eq 2 then $
          symm_triangles,xx,yy,triangles,$
                         triangconj,ntriangles,rectangles,nrectan

      xreg=dblarr(nxreg(0),nxreg(1),2)
      dx=(xreglimits(2)-xreglimits(0))/(nxreg(0)-1)
      dy=(xreglimits(3)-xreglimits(1))/(nxreg(1)-1)
      for i=0,nxreg(1)-1 do xreg(*,i,0)=dx*indgen(nxreg(0))+xreglimits(0)
      for i=0,nxreg(0)-1 do xreg(i,*,1)=dy*indgen(nxreg(1))+xreglimits(1)

   endif
   if n_elements(wregpad) ne nw then begin
      wregpad=dblarr(nw)
      for iw=0,nw-1 do begin
         wmax=max(w(*,*,iw))
         wmin=min(w(*,*,iw))
         if wmax*wmin lt 0 then wregpad(iw)=0 $
         else                   wregpad(iw)=wmin-0.1*(wmax-wmin)
      endfor
   endif

   wreg=dblarr(nxreg(0),nxreg(1),nw)

   case 1 of

   symmtri eq 3: for iw=0,nw-1 do $
      wreg(*,*,iw)=grid_data(xx,yy,reform(w(*,*,iw)),nxreg,xreglimits,$
                             triangles,wregpad(iw))

   symmtri eq 0 or (symmtri lt 3 and nrectan eq 0): for iw=0,nw-1 do $
      wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
           [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

   symmtri eq 1 and nrectan gt 0: $
      fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
           triangles,ntriangles,rectangles

   symmtri eq 2 and nrectan gt 0: $
      average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
           triangles,triangconj

   endcase

   err=check_math(0,0)
   ;Floating underflow is not a real error, the message is suppressed
   if err ne 32 and err ne 0 then print,'Math error in regulargrid:',err

end

;==============================================================================
pro symm_triangles,xx,yy,triangles,$
         triangconj,ntriangles,rectangles,nrectan

   ntriangles=n_elements(triangles(0,*))
   print,'Triangulation includes ',ntriangles, '  triangles'
   print,'Checking triangulation ...'

   npoints=n_elements(xx)

   dist=dblarr(npoints-1)
   for i=0L,npoints-2 do $
      dist(i)=(xx(i+1)-xx(i))^2+(yy(i+1)-yy(i))^2
   dist2=min(dist)
   rectangles=lonarr(3,ntriangles)
   ;Structure of the rectangles array:
   ;If(rectangles(0,i)=1 then the Ith triangle from the triangles array
   ;is rectangular one
   tmp1=lonarr(3) & nrec_tri=0
   for i=0L,ntriangles-1 do begin
      if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
         (xx(triangles(1,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(1,i)))*$
         (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
	 rectangles(0,i)=1
         tmp1(0)=triangles(1,i)
         if xx(triangles(0,i)) lt xx(triangles(2,i)) then begin
            tmp1(1)=triangles(0,i)
            tmp1(2)=triangles(2,i)
         endif else begin
            tmp1(1)=triangles(2,i)
            tmp1(2)=triangles(0,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif

      if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
         (xx(triangles(0,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(1,i)))*$
         (yy(triangles(0,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
         rectangles(0,i)=1
         tmp1(0)=triangles(0,i)
         if xx(triangles(1,i)) lt xx(triangles(2,i)) then begin
            tmp1(1)=triangles(1,i)
            tmp1(2)=triangles(2,i)
         endif else begin
            tmp1(1)=triangles(2,i)
            tmp1(2)=triangles(1,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif

      if abs((xx(triangles(0,i))-xx(triangles(2,i)))*$
         (xx(triangles(1,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(2,i)))*$
         (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
         rectangles(0,i)=1
         tmp1(0)=triangles(2,i)
         if xx(triangles(0,i)) lt xx(triangles(1,i)) then begin
            tmp1(1)=triangles(0,i)
            tmp1(2)=triangles(1,i)
         endif else begin
            tmp1(1)=triangles(1,i)
            tmp1(2)=triangles(0,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif
   endfor
   ;Rectangles(1,i) is not equal to zero if the ith rectangular triandgle
   ;has a common long side with the jth rectangular triangle. In this case
   ;rectangles(2,i)=j
   nrectan=0L
   for i=0L,ntriangles-1 do begin
     if rectangles(0,i) gt 0 then begin
        nrec_tri=nrec_tri+1
        if rectangles(1,i) eq 0 then begin
        for j=i+1L,ntriangles-1 do begin
           if rectangles(0,j) gt 0 then $
           if triangles(1,i) eq triangles(1,j) then $
           if triangles(2,i) eq triangles(2,j) then begin
              nrectan=nrectan+1
              rectangles(1,i)=1
              rectangles(2,i)=j
              rectangles(1,j)=1
              rectangles(2,j)=i
              goto,out
           endif
        endfor
        out:
        endif
     endif
   endfor

   if nrectan ne 0  then begin
      print,'Among    ',nrec_tri, '  rectangular triangles'
      print,'there are',nrectan, '   pairs which have common circumcircle'
      tmp2=lonarr(4)
      ndiag1=0
      ndiag2=0
      triangconj=lonarr(3,ntriangles)
      for i=0L,ntriangles-1 do begin
         if rectangles(1,i) eq 1 then begin
            if rectangles(2,i) gt i then begin
               for j=0,2 do tmp2(j)=triangles(j,i)
               tmp2(3)=triangles(0,rectangles(2,i))
               if yy(tmp2(1)) lt yy(tmp2(2)) then ndiag1=ndiag1+1 else $
                  ndiag2=ndiag2+1
               triangconj(0,i)=tmp2(1)
               triangconj(1,i)=tmp2(0)
               triangconj(2,i)=tmp2(3)
               triangconj(0,rectangles(2,i))=tmp2(2)
               triangconj(1,rectangles(2,i))=tmp2(0)
               triangconj(2,rectangles(2,i))=tmp2(3)
            endif
         endif else for j=0,2 do triangconj(j,i)=triangles(j,i)
      endfor
      print,' Among them ',ndiag1, ' are formed by the triangles,'
      print,' having the common side which is oriented as /////'
      print,' and ',ndiag2, ' have the triangles common side,'
      print,' oriented as \\\\\'
      print,' Calculating the conjugated triangulation ...'
   endif

end

;==============================================================================
pro fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
        triangles,ntriangles,rectangles

   tmp2=lonarr(4)

   for iw=0,nw-1 do begin
      data=reform(w(*,*,iw))
      print,'Calculating the fitting triangulization for iw=',iw
      for i=0L,ntriangles-1 do begin
         if rectangles(1,i) eq 1 then begin
            if rectangles(2,i) gt i then begin
               tmp2(0:2)=triangles(0:2,i)
               tmp2(3)  =triangles(0,rectangles(2,i))
               if abs(data(tmp2(0))-data(tmp2(3))) lt $
                  abs(data(tmp2(1))-data(tmp2(2))) then begin
                  triangles(0,i)=tmp2(1)
                  triangles(1,i)=tmp2(0)
                  triangles(2,i)=tmp2(3)
                  triangles(0,rectangles(2,i))=tmp2(2)
                  triangles(1,rectangles(2,i))=tmp2(0)
                  triangles(2,rectangles(2,i))=tmp2(3)
               endif
            endif
         endif
      endfor
      wreg(*,*,iw)=trigrid(xx,yy,data,triangles, $
         [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))
   endfor
   print,'Using fitted triangulation'
end

;==============================================================================
pro average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
        triangles,triangconj

   wconj=dblarr(nxreg(0),nxreg(1))

   for iw=0,nw-1 do begin
      ; Calculate wreg with original triangulation
      wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
            [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

      ; Calculate wconj with conjugated triangulation
      wconj       =trigrid(xx,yy,w(*,*,iw),triangconj, $
            [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

      wreg(*,*,iw) = 0.5*(wreg(*,*,iw) + wconj)
   endfor
   print,'Using averaged conjugated triangulation'

end
;===========================================================================
pro polargrid,nvector,vectors,x,w,xreg,wreg
;
;  Transform 2D or 3D grid and vector variables 
;  from Cartesian to cylindrical components

  on_error,2

  siz = size(x)
  ndim = siz(0)-1
  case ndim of
      2: polargrid2,nvector,vectors,x,w,xreg,wreg
      3: polargrid3,nvector,vectors,x,w,xreg,wreg
      else: begin
          print,'polargid works for 2D and 3D arrays only'
          retall
      end
  endcase
end
;===========================================================================
pro polargrid2,nvector,vectors,x,w,xreg,wreg
;
;  Transform 2D grid and vector variables 
;  from x and y to radial and phi components

  on_error,2

  xreg=x
  xreg(*,*,0)=sqrt(x(*,*,0)^2+x(*,*,1)^2)
  xreg(*,*,1)=atan(x(*,*,1),x(*,*,0))
  phi=xreg(*,*,1)
  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)+w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)= -w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor

  ;Remove 2*pi jumps from phi
  pi2=8*atan(1) & sz=size(phi) & nx2=sz(2)
  for ix2=1,nx2-1 do while phi(1,ix2-1) gt phi(1,ix2) do $
     phi(*,ix2)=phi(*,ix2)+pi2

  xreg(*,*,1)=phi
end

;===========================================================================
pro polargrid3,nvector,vectors,x,w,xreg,wreg
;
;    Transform 3D vector variables from x,y,z to radial,phi,z components

  on_error,2

  xreg=x
  xreg(*,*,*,0)=sqrt(x(*,*,*,0)^2+x(*,*,*,1)^2)
  xreg(*,*,*,1)=atan(x(*,*,*,1),x(*,*,*,0))
  phi=xreg(*,*,*,1)
  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,*,ivx)=  w(*,*,*,ivx)*cos(phi)+w(*,*,*,ivy)*sin(phi)
     wreg(*,*,*,ivy)= -w(*,*,*,ivx)*sin(phi)+w(*,*,*,ivy)*cos(phi)
  endfor

  ;Remove 2*pi jumps from phi
  pi2=8*atan(1) & sz=size(phi) & nx2=sz(2)
  for ix2=1,nx2-1 do while phi(1,ix2-1) gt phi(1,ix2) do $
     phi(*,ix2)=phi(*,ix2)+pi2

  xreg(*,*,*,1)=phi
end

;===========================================================================
pro spheregrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x,y,z to radial,phi,z components

  on_error,2

  xreg=x
  xreg(*,*,*,0)=sqrt(x(*,*,*,0)^2+x(*,*,*,1)^2+x(*,*,*,2)^2)
  xreg(*,*,*,2)=-atan(x(*,*,*,2),x(*,*,*,0))
  xreg(*,*,*,1)=atan(x(*,*,*,1),sqrt(x(*,*,*,0)^2+x(*,*,*,2)^2))
  phi=xreg(*,*,*,2)
  theta=xreg(*,*,*,1)
  wreg=w
  sinphi=sin(phi)
  cosphi=cos(phi)
  sintheta=sin(theta)
  costheta=cos(theta)
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     ivz=ivy+1
     wreg(*,*,*,ivx)=(w(*,*,*,ivx)*cosphi-w(*,*,*,ivz)*sinphi)*costheta $
                     +w(*,*,*,ivy)*sintheta
     wreg(*,*,*,ivz)=-w(*,*,*,ivx)*sinphi-w(*,*,*,ivz)*cosphi
     wreg(*,*,*,ivy)=(-w(*,*,*,ivx)*cosphi+w(*,*,*,ivz)*sinphi)*sintheta $
                     +w(*,*,*,ivy)*costheta
  endfor

  ;Remove 2*pi jumps from phi
  pi=4*atan(1) & pi2=2*pi & sz=size(phi) & nx2=sz(2) & nx3=sz(3)
  for ix3=1,nx3-1 do while phi(1,1,ix3-1) gt phi(1,1,ix3) do $
     phi(*,*,ix3)=phi(*,*,ix3)+pi2

  ;Remove turn over from theta
  for ix2=1,nx2-1 do $
  if theta(1,ix2-1,1) ge theta(1,ix2,1) then begin
     if theta(1,ix2,1) lt 0 then $
          theta(*,ix2-1,*)=-pi-theta(*,ix2-1,*) $
     else $
          theta(*,ix2,*)=pi-theta(*,ix2,*)
  endif

  xreg(*,*,*,2)=phi
  xreg(*,*,*,1)=theta
end

;===========================================================================
pro unpolargrid,nvector,vectors,x,w,xreg,wreg
;
; Transform 2D and 3D vector variables cylindrical to Cartesian components

  on_error,2

  siz = size(x)
  ndim = siz(0)-1
  case ndim of
      2: unpolargrid2,nvector,vectors,x,w,xreg,wreg
      3: unpolargrid3,nvector,vectors,x,w,xreg,wreg
      else: begin
          print,'unpolargid works for 2D and 3D arrays only'
          retall
      end
  endcase

end
;===========================================================================
pro unpolargrid2,nvector,vectors,x,w,xreg,wreg
;
;    Transform 2D grid (and vector variables)
;    from radial and phi to x and y components

  xreg=x
  phi=x(*,*,1)

  phimax = max(abs(phi))

  ; If phi is in local time or degrees, change it to radians
  if      abs(phimax - 24.0) lt 0.1 then phi=phi*!pi/12 $
  else if phimax gt 20.0            then phi=phi*!pi/180

  xreg(*,*,0)=x(*,*,0)*cos(phi)
  xreg(*,*,1)=x(*,*,0)*sin(phi)

  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)-w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)=  w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor
end

;===========================================================================
pro unpolargrid3,nvector,vectors,x,w,xreg,wreg
;
;  Transform 3D grid (and vector variables)
;  from radial, phi,z to x, y, z components

  on_error,2

  xreg=x
  phi=x(*,*,*,1)

  ; If phi is in degrees, change it to radians
  if max(abs(phi)) gt 20. then phi=phi*!pi/180

  xreg(*,*,*,0)=x(*,*,*,0)*cos(phi)
  xreg(*,*,*,1)=x(*,*,*,0)*sin(phi)

  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,*,ivx)=  w(*,*,*,ivx)*cos(phi)-w(*,*,*,ivy)*sin(phi)
     wreg(*,*,*,ivy)=  w(*,*,*,ivx)*sin(phi)+w(*,*,*,ivy)*cos(phi)
  endfor
end

;===========================================================================
pro getaxes,ndim,x,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
;===========================================================================
on_error,2
case ndim of
  1: xx=x
  2: begin
        xx=x(*,*,0)
        yy=x(*,*,1)
     end
  3: begin
       xx=x(*,*,*,0)
       yy=x(*,*,*,1)
       zz=x(*,*,*,2)
     end
endcase

if keyword_set(cut0) then begin
   xx=xx(cut0)
   if ndim gt 1 then yy=yy(cut0)
   if ndim gt 2 then zz=zz(cut0)
endif

!x.title=strupcase(variables(0))
if plotdim gt 1 then !y.title=strupcase(variables(1))
if plotdim gt 2 then !z.title=strupcase(variables(2))

if ndim eq 3 and plotdim eq 2 then begin
   siz=size(cut)
   case 1 of
     siz(0) eq 2: rSlice=zz(0)
     siz(1) eq 1: rSlice=xx(0)
     siz(2) eq 1: rSlice=yy(0)
   endcase
   print,'Normal coordinate of 2D slice:',rSlice
endif else        rSlice=0.0

; Cut with fixed X value?
siz=size(cut)
; in 2D
if siz(0) eq 2 and siz(1) eq 1 then begin
   xx=yy
   !x.title=variables(1)
endif
; in 3D
if siz(0) eq 3 then begin
   case 1 of
   plotdim eq 1: begin
      xx=zz
      !x.title=variables(2)
   end
   siz(1) eq 1: begin
         xx=yy
         yy=zz
         !x.title=variables(1)
         !y.title=variables(2)
   end
   siz(2) eq 1: begin
      yy=zz
      !y.title=variables(2)
   end
   else: print,'internal error in getaxes'
   endcase
endif

end

;===========================================================================
pro getfunc,f,f1,f2,func1,func2,x,w,physics,eqpar,variables,cut,rcut
;===========================================================================
on_error,2

f1=funcdef(x,w,func1,physics,eqpar,variables,rcut)

if keyword_set(cut) then f1=f1(cut)

if func2 eq '' then f=f1 else begin

   f2=funcdef(x,w,func2,physics,eqpar,variables,rcut)

   if keyword_set(cut) then f2=f2(cut)

   ; Calculate f=sqrt(f1^2+f2^2)
   f=sqrt(f1^2+f2^2)
endelse

end

;===========================================================================
pro plot_func,x,w,xreg,wreg,usereg,ndim,physics,eqpar,rBody,$
  variables,wnames,axistype,plotmodes,plottitles,$
  ax,az,contourlevel,linestyle,$
  velvector,velspeed,velseed,velpos,velx,vely,veltri,$
  cut,cut0,rcut,plotdim,$
  nfunc,multix,multiy,fixaspect,plotix,plotiy,funcs,funcs1,funcs2,fmin,fmax,f
;===========================================================================
   on_error,2

   ; Get grid dimensions and set irr=1 if it is an irregular grid

   if keyword_set(cut) then siz = size(cut)  $
   else if usereg then      siz = size(xreg) $
   else                     siz = size(x)
   nx=siz(1)
   if plotdim eq 1 then begin
      ny=1
      irr=0
   endif else begin
      ny=siz(2)
      irr=ny eq 1
   endelse

   if irr and axistype eq 'cells' then begin
      print,'Irregular grid, axistype must be set to coord'
      axistype='coord'
   endif

   if axistype eq 'coord' then begin
      if usereg then $
         getaxes,ndim,xreg,xx,yy,zz,cut,cut0,rSlice,plotdim,variables $
      else $
         getaxes,ndim,x   ,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
   endif

   ; Calculate plot spacing from number of plots per page (ppp) and charsize
   if !p.charsize eq 0.0 then !p.charsize=1.0
   ppp   = multix*multiy
   space = max([float(!d.y_ch_size)/float(!d.y_size),$
                float(!d.x_ch_size)/float(!d.x_size)])*3.0*!p.charsize
   set_space, ppp, space, sizes, nx = multix, ny = multiy

   ; Store x and y titles and tick names
   xtitle    = !x.title
   ytitle    = !y.title
   xtickname = !x.tickname
   ytickname = !y.tickname

   for ifunc=0,nfunc-1 do begin

      plotmod=plotmodes(ifunc)

      ; stream2 --> stream
      i=strpos(plotmod,'stream2')
      if i ge 0 then plotmod=strmid(plotmod,0,i+6)+strmid(plotmod,i+7)

      ; contour --> cont
      i=strpos(plotmod,'contour')
      if i ge 0 then plotmod=strmid(plotmod,0,i+4)+strmid(plotmod,i+7)

      i=strpos(plotmod,'white')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+5)
          white=1
      endif else white=0

      i=strpos(plotmod,'grid')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showgrid=1
      endif else showgrid=0

      i=strpos(plotmod,'irr')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
          irr=1
      endif

      i=strpos(plotmod,'mesh')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showgrid=1
          if irr then showmesh=0 else showmesh=1
      endif else showmesh=0

      i=strpos(plotmod,'body')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showbody=1
      endif else showbody=0

      i=strpos(plotmod,'fill')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          fill=1
      endif else fill=0

      i=strpos(plotmod,'bar')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
          showbar=1
          if strpos(plotmod,'stream') ge 0 then showbar=0
          fill=1
      endif else showbar=0

      i=strpos(plotmod,'label')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+5)
          label=1
      endif else label=0

      i=strpos(plotmod,'over')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          !p.multi(0) = !p.multi(0)+1
          if !p.multi(0) ge !p.multi(1)*!p.multi(2) then !p.multi(0)=0
      endif

      i=strpos(plotmod,'log')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
          logarithm=1
      endif else logarithm=0

      !p.title=plottitles(ifunc)
      if !p.title eq 'default' then !p.title=funcs(ifunc)

      ; Calculate the next p.multi(0) explicitly
      if !p.multi(0) gt 0 then multi0=!p.multi(0)-1 $
      else multi0=!p.multi(1)*!p.multi(2)-1

      ; Calculate subplot position indices
      if !p.multi(4) then begin
         ; columnwise
         plotix=multix-1-multi0/multiy
         plotiy=multi0 mod multiy
      endif else begin
         ; rowwise
         plotix=multix-1-(multi0 mod multix)
         plotiy=multi0/multix
      endelse

      if plotmod ne 'shade' and plotmod ne 'surface' then begin

        ; obtain position for flat plotmodes
        set_position, sizes, plotix, multiy-1-plotiy, pos, /rect

        ; shrink in X direction for (overplotting) a colorbar
        if strpos(plotmodes(ifunc),'bar') ge 0 $
          or (strpos(plotmodes(ifunc),'over') ge 0 and $
              strpos(plotmodes(ifunc-1>0),'bar') ge 0) then $
          pos(2) = pos(2) - (pos(2) - pos(0))*0.15

        ; shrink in X direction for the Y axis of plot
        if plotmod eq 'plot' and multix gt 1 then $
          pos(0) = pos(0) + (pos(2) - pos(0))*0.15

        if keyword_set(fixaspect) and plotmod ne 'plot' then begin

	  if plotmod eq 'polar' then $
            aspectx=1 $
          else begin
            if !x.range(1) ne !x.range(0) then    $
               width=abs(!x.range(1)-!x.range(0)) $
            else if axistype eq 'coord' then      $
               width=  max(xx) - min(xx)          $
            else                                  $
               width=  nx-1.0

            if !y.range(1) ne !y.range(0) then    $
               height=abs(!y.range(1)-!y.range(0))$
            else if axistype eq 'coord' then      $
               height= max(yy) - min(yy)          $
            else                                  $
               height= ny-1.0

            aspectx = width/height
          endelse

          aspectpos = (pos(2)-pos(0))/(pos(3)-pos(1)) $
                    *float(!d.x_size)/float(!d.y_size)

          aspectratio = aspectpos/aspectx

          ;print,'aspectx,pos,ratio=',aspectx,aspectpos,aspectratio

          if aspectratio gt 1 then begin
             posmid=(pos(2)+pos(0))/2.
             posdif=(pos(2)-pos(0))/2.
             pos(0)=posmid - posdif/aspectratio
             pos(2)=posmid + posdif/aspectratio
          endif else begin
             posmid=(pos(3)+pos(1))/2.
             posdif=(pos(3)-pos(1))/2.
             pos(1)=posmid - posdif*aspectratio
             pos(3)=posmid + posdif*aspectratio
          endelse
        endif

        ; Omit X axis if unneeded
        if (plotiy gt 0) then begin
          !x.tickname = strarr(60)+' '
          !x.title = ' '
        endif

        ; Omit Y axis if unneeded
        if (plotix gt 0 and plotmod ne 'plot') then begin
          !y.tickname = strarr(60)+' '
          !y.title = ' '
        endif

        !p.position = pos

      endif

      if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             xreg,wreg,physics,eqpar,variables,cut0,rcut $
      else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             x,  w,   physics,eqpar,variables,cut0,rcut

      f_min=fmin(ifunc)
      f_max=fmax(ifunc)

      if logarithm and f_min gt 0 and min(f) gt 0 then begin
          f     = alog10(f)
          f_min = alog10(f_min)
          f_max = alog10(f_max)
          if plottitles(ifunc) eq 'default' then !p.title = 'log '+!p.title
      endif

      if f_max eq f_min then begin
         f_max=f_max+1
         f_min=f_min-1
      endif

      if plotmod eq 'plot' then $
         if nfunc gt ppp                then lstyle=ifunc/ppp $
         else if keyword_set(linestyle) then lstyle=linestyle $
         else                                lstyle=!p.linestyle

      ; Skip minimum ad maximum levels
      if plotmod eq 'cont' or plotmod eq 'polar' then $
         levels=(findgen(contourlevel+2)-1)/(contourlevel-1) $
                *(f_max-f_min)+f_min

      ; figure out the units of angle in the second coordinate
      if plotmod eq 'polar' then begin
          if max(yy)-min(yy) gt 300 then $
            angleunit = !pi/180 $   ; degrees
          else if max(yy)-min(yy) gt 20 then $
            angleunit = !pi/12 $    ; local time
          else $
            angleunit = 1.0         ; radians
      endif

      if plotmod eq 'tv' then begin
         ; Calculate plotting position and size

         tvplotx=pos(0)*!d.x_size
         tvploty=pos(1)*!d.y_size
         tvsizex=(pos(2)-pos(0))*!d.x_size
         tvsizey=(pos(3)-pos(1))*!d.y_size
         ; recalculate f for tv mode
         if !d.name eq 'PS' then tvf=congrid(f,200,200) $
         else                    tvf=congrid(f,tvsizex,tvsizey)

         tvf=bytscl(tvf,MIN=f_min,MAX=f_max,TOP=!D.TABLE_SIZE-3)+1
      endif

      if showbar then plot_color_bar, $
        [pos(2)+0.005, pos(1), pos(2)+0.025, pos(3)], [f_min,f_max]

      case axistype of
      'cells': case plotmod of
         'cont': contour,f>f_min,LEVELS=levels,$
                 FILL=fill,FOLLOW=label,XSTYLE=1,YSTYLE=1,/NOERASE
         'plot'     :plot,f,YRANGE=[f_min,f_max],XSTYLE=18,ystyle=18, $
                                                 LINE=lstyle,/NOERASE
         'shade'    :begin
                        shade_surf,f>f_min,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                        if showgrid then $
                           surface,f>f_min,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                     end
         'surface'  :surface,f>f_min,ZRANGE=[f_min,f_max],$
                        XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :begin
                        tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
                        contour,f,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
                     end
         'vel'      :vector,f1,f2,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'vector'   :vector,f1,f2,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,$
                        /NOERASE,WHITE=white
         'stream'   :streamline,f1,f2,NVECS=velvector,SEED=velseed,X0=velpos,$
                        /NOERASE,WHITE=white
         'velovect' :velovect,f1,f2,/NOERASE
         'ovelovect':velovect,f1,f2,/NOERASE,$
            XRANGE=[0,n_elements(f1(*,0))-1],YRANGE=[0,n_elements(f1(0,*))-1]
         endcase
      'coord': case plotmod of
         'cont'     :if irr then begin
                       if not keyword_set(tri) then $
                         triangulate,float(xx),float(yy),tri
                       contour,f>f_min,xx,yy,$
                          FOLLOW=label, FILL=fill, TRIANGULATION=tri, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
                    endif else $
                       contour,f>f_min,xx,yy,$
                          FOLLOW=label, FILL=fill, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
	 'polar'    :polar_contour,f>f_min,yy*angleunit,xx,$
                          FOLLOW=label, FILL=fill, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
         'plot'     :plot,xx,f,YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,$
                          LINE=lstyle,/NOERASE
         'shade'    :if irr then begin
                        shade_surf_irr,f>f_min,xx,yy,AX=ax,AZ=az
                        shade_surf,f>f_min,xx,yy,AX=ax,AZ=az,/NODATA,/NOERASE
                     endif else begin
                        shade_surf,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                        if showgrid then $
                           surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                     endelse
         'surface'  :surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                        XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :begin
                       tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
                       contour,f,xx,yy,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
                     end
         'vel'      :vector,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                        TRIANGLES=veltri,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'vector'   :vector,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                        TRIANGLES=veltri,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,$
                        /NOERASE, WHITE=white
         'stream'   :streamline,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                        TRIANGLES=veltri,NVECS=velvector,$
                        SEED=velseed,X0=velpos,$
                        /NOERASE, WHITE=white
         'velovect' :velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE
         'ovelovect':velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE,$
                        XRANGE=[min(xx),max(xx)],YRANGE=[min(yy),max(yy)]
         endcase
      else:print,'Unknown axistype:',axistype
      endcase

      if showbody and axistype eq 'coord' then $
      if rBody gt abs(rSlice) then begin
         theta = findgen(37)*!pi*2.0/36.0
         rBodySlice=sqrt(rBody^2-rSlice^2)
         polyfill, rBodySlice*cos(theta), rBodySlice*sin(theta),color = 0, $
           noclip=0
         ; redraw box in case the body is at the edge
         plot,xx,yy,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
      endif

      if showgrid and plotdim eq 2 and plotmod ne 'surface'    $
                                   and plotmod ne 'shade' then begin
          if(plotmod eq 'polar')then                                       $
            plotgrid,xx,yy*angleunit,lines=showmesh,xstyle=1,ystyle=1,/polar $
          else if keyword_set(cut) then                                    $
            plotgrid,xx,yy,lines=showmesh,xstyle=1,ystyle=1                $
          else begin
              if !x.range[0] ne !x.range[1] then xrange=!x.range else $
                xrange=[min(xx),max(xx)]
              if !y.range[0] ne !y.range[1] then yrange=!y.range else $
                yrange=[min(yy),max(yy)]
              plotgrid,x,lines=showmesh,xstyle=1,ystyle=1,$
                xrange=xrange,yrange=yrange
          endelse
      endif

      !p.multi(0) = multi0
      !p.position = 0
      !x.title    = xtitle
      !x.tickname = xtickname
      !y.title    = ytitle
      !y.tickname = ytickname
   endfor

   !p.position = 0

end
;===========================================================================
pro putbottom,multix,multiy,ix,iy,info,nx,it,time

on_error,2

t    = round(time)
sec  = t mod 60
t    = t/60
min  = t mod 60
t    = t/60
hour = t mod 24
t    = t/24
day  = t mod 365
year = t/365

if year gt 0 then $
  stime = 'time='+string(year,day,hour,format='(i3,"y",i3.3,"d",i2.2,"h")') $
else if day gt 0 then $
  stime = 'time='+string(day,hour,min,format='(i4,"d",i2.2,"h",i2.2,"m")') $
else if hour gt 0 then $
  stime = 'time='+string(hour,min,sec,format='(i4,"h",i2.2,"m",i2.2,"s")') $
else $
  stime = 'time='+string(time, format='(g12.5)')

snx = string(nx, format='(i6,2(",",i4))')
sit = string(it, format='(i8)')

siz = size(info)

if siz(1) eq 7 then begin
  if not execute("s = "+info) then begin
      print,'Error in putbottom: cannot evaluate info=',info
      return
  endif
endif else begin
    if info lt 1 then return
    s = ''
    if info gt 2 then s = 'nx=' + snx + ' '
    if info gt 1 then s = s + 'it=' + sit + ', '
    s = s + stime
endelse

xyouts, 5+(ix*!d.x_size)/multix, 8+(iy*!d.y_size)/multiy, /DEV, s

;;xyouts,5,3000,/DEV,info

;info='time ='+string(time/(17.24*3600),format='(f6.2)')+' Uday'
;xyouts,3000+(0.9*ix*!d.x_size)/multix,3000+(0.82*(iy+1)*!d.y_size)/multiy,/DEV,info,charsize=1.2

;info='time ='+string(time/60,format='(i4)')+' min'
;xyouts,!p.position[0]-0.1*(!p.position(2)-!p.position(0)),$
;       !p.position[1]-0.1*(!p.position(3)-!p.position(1)),$
;       info,/NORM,charsize=0.8

end
;===========================================================================
pro putheader,multix,multiy,ix,iy,ninfo,headline,nx

on_error,2

if ninfo lt 1 then return
info=strtrim(headline,2)
if ninfo gt 1 then info=info+' (nx='+string(nx,format='(i6,2(i4))')+')'
xyouts,5+(ix*!d.x_size)/multix,-12+((iy+1)*!d.y_size)/multiy,/DEV,info

end
;===========================================================================
function diff1,a,x
;
; Take derivative of "a" with respect to "x" (if present)
; using 2nd order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
ndim=siz(0)
if ndim ne 1 then begin
   print,'Function diff1 is intended for 1D arrays only'
   retall
endif

n  = n_elements(a)
nx = n_elements(x)

if nx ne 0 and nx ne n then begin
   print,'Error in diff1, arrays sizes differ: nx, n=', nx, n
   retall
endif

dadx = a

if nx eq 0 then dadx(1:n-2) = (a(2:n-1) - a(0:n-3))/2 $
else            dadx(1:n-2) = (a(2:n-1) - a(0:n-3))/(x(2:n-1) - x(0:n-3))

; fill in boundaries
dadx(0)   = dadx(1)
dadx(n-1) = dadx(n-2)

return,dadx
end
;===========================================================================
function laplace,a,x,y
;
; Take Laplace of "a" in 1 or 2D with respect to "x" (and "y") if present.
;
;===========================================================================
on_error,2

siz=size(a)
ndim=siz(0)
if ndim eq 2 then return,laplace2(a,x,y)
if ndim ne 1 then begin
   print,'Function laplace is intended for 1D and 2D arrays only'
   retall
endif

n  = n_elements(a)
nx = n_elements(x)

if nx ne 0 and nx ne n then begin
   print,'Error in diff1, arrays sizes differ: nx, n=', nx, n
   retall
endif

d2adx2 = a

if nx eq 0 then d2adx2(1:n-2) = a(2:n-1) - 2*a(1:n-2) + a(0:n-3) $
else            d2adx2(1:n-2) = $
  ( (a(2:n-1) - a(1:n-2))/(x(2:n-1) - x(1:n-2)) $
  - (a(1:n-2) - a(0:n-3))/(x(1:n-2) - x(0:n-3)) ) $
  / (0.5*(x(2:n-1) - x(0:n-3)))

; fill in boundaries
d2adx2(0)   = d2adx2(1)
d2adx2(n-1) = d2adx2(n-2)

return,d2adx2
end
;===========================================================================
function laplace2,a,x,y
;
; Take Laplace of "a" in 2D with respect to "x" and "y" (if present).
;
;===========================================================================
on_error,2

siz=size(a)
ndim=siz(0)

n1 = siz(1)
n2 = siz(2)
n  = n1*n2
nx = n_elements(x)
ny = n_elements(y)

if (nx ne 0 or ny ne 0) and (nx ne n or ny ne n) then begin
   print,'Error in diff1, arrays sizes differ: n, nx, ny=', n, nx, ny
   retall
endif

d2adx2 = a

if nx eq 0 then d2adx2(1:n1-2,1:n2-2) = $
  a(2:n1-1,1:n2-2) - 2*a(1:n1-2,1:n2-2) + a(0:n1-3,1:n2-2) + $
  a(1:n1-2,2:n2-1) - 2*a(1:n1-2,1:n2-2) + a(1:n1-2,0:n2-3)   $
else            d2adx2(1:n1-2,1:n2-2) = $
  ( (a(2:n1-1,1:n2-2) - a(1:n1-2,1:n2-2))/ $
    (x(2:n1-1,1:n2-2) - x(1:n1-2,1:n2-2))- $
    (a(1:n1-2,1:n2-2) - a(0:n1-3,1:n2-2))/ $
    (x(1:n1-2,1:n2-2) - x(0:n1-3,1:n2-2)) $
  ) / (0.5*(x(2:n1-1,1:n2-2) - x(0:n1-3,1:n2-2))) + $
  ( (a(1:n1-2,2:n2-1) - a(1:n1-2,1:n2-1))/ $
    (y(1:n1-2,2:n2-1) - y(1:n1-2,1:n2-1))- $
    (a(1:n1-2,1:n2-2) - a(1:n1-2,0:n2-3))/ $
    (y(1:n1-2,1:n2-2) - y(1:n1-2,0:n2-3)) $
  ) / (0.5*(y(1:n1-2,2:n2-1) - y(1:n1-2,0:n2-3)))

; fill in boundaries
d2adx2(0   ,1:n2-1) = d2adx2(1   ,1:n2-1)
d2adx2(n1-1,1:n2-1) = d2adx2(n1-2,1:n2-1)
d2adx2(*,0)         = d2adx2(*,1)
d2adx2(*,n2-1)      = d2adx2(*,n2-2)

return,d2adx2
end
;===========================================================================
function diff2,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 2nd order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
ndim=siz(0)
if ndim ne 2 and ndim ne 3 then begin
   print,'Function diff2 is intended for 2D and 3D arrays only'
   retall
endif

if direction lt 1 or direction gt ndim then begin
   print,'Direction=',direction,' should be between 1 and ndim=',ndim,'!'
   retall
endif

n1=siz(1)
n2=siz(2)
if ndim eq 3 then n3=siz(3)

if direction eq 1 then begin
   ind1=indgen(n1)
   jnd1=ind1+1
   jnd1(n1-1)=n1
   hnd1=ind1-1
   hnd1(0)=0
   if ndim eq 2 then $
     dadx=(a(jnd1,*)-a(hnd1,*))/(x(jnd1,*)-x(hnd1,*)) $
   else $
     dadx=(a(jnd1,*,*)-a(hnd1,*,*))/(x(jnd1,*,*)-x(hnd1,*,*))
endif
if direction eq 2 then begin
   ind2=indgen(n2)
   jnd2=ind2+1
   jnd2(n2-1)=n2
   hnd2=ind2-1
   hnd2(0)=0
   if ndim eq 2 then $
     dadx=(a(*,jnd2)-a(*,hnd2))/(x(*,jnd2)-x(*,hnd2)) $
   else $
     dadx=(a(*,jnd2,*)-a(*,hnd2,*))/(x(*,jnd2,*)-x(*,hnd2,*))
endif
if direction eq 3 then begin
   ind3=indgen(n3)
   jnd3=ind3+1
   jnd3(n2-1)=n3
   hnd3=ind3-1
   hnd3(0)=0
   dadx=(a(*,*,jnd3)-a(*,*,hnd3))/(x(*,*,jnd3)-x(*,*,hnd3))
endif

return,dadx

end

;===========================================================================
function diff4,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 4th order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff4 is intended for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

dadx=a

if direction eq 1 then begin
   if n1 lt 5 then begin
      print,'Cannot take 4th order X gradient of grid with less than 5 columns'
      retall
   endif
   dadx(2:n1-3,*)=(a(4:n1-1,*)-8*a(3:n1-2,*)+8*a(1:n1-4,*)-a(0:n1-5,*)) $
                 /(x(3:n1-2,*)-x(1:n1-4,*))/6;
   dadx(0,*)   =dadx(2,*)
   dadx(1,*)   =dadx(2,*)
   dadx(n1-2,*)=dadx(n1-3,*)
   dadx(n1-1,*)=dadx(n1-3,*)
endif
if direction eq 2 then begin
   if n2 lt 5 then begin
      print,'Cannot take 4th order Y gradient of grid with less than 5 rows'
      retall
   endif
   dadx(*,2:n2-3)=(a(*,4:n2-1)-8*a(*,3:n2-2)+8*a(*,1:n2-4)-a(*,0:n2-5)) $
                 /(x(*,3:n2-2)-x(*,1:n2-4))/6;
   dadx(*,0)   =dadx(*,2)
   dadx(*,1)   =dadx(*,2)
   dadx(*,n2-2)=dadx(*,n2-3)
   dadx(*,n2-1)=dadx(*,n2-3)
endif

return,dadx

end

;===========================================================================
function diff3,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using IDL's 1D deriv() function
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff3 is intended for 2D arrays only'
   retall
endif

dadx=a

if direction eq 1 then for i2=0,siz(2)-1 do dadx(*,i2)=deriv(x(*,i2),a(*,i2))
if direction eq 2 then for i1=0,siz(1)-1 do dadx(i1,*)=deriv(x(i1,*),a(i1,*))

return,dadx

end

;===========================================================================
function minmod,a,b
;
; Calculate minmod limited slope of a and b slopes

on_error,2

; get sign of a
if a gt 0 then s=1 else s=-1

; calculate limited slope
c = s*max([0,min([abs(a),s*b])])

return,c
end
;==========================================================================
function symmdiff,direction,a,x,y,report=report,anti=anti
;
;find the symmetric for irregular grid
;
;=========================================================================

; If x is not present, the grid is taken to be Cartesian
if n_elements(x) eq 0 then return, symmdiffreg(direction,a)

if not keyword_set(report) then report = 0

diff=a*0

e = exp(1.0d0)
sorting = x + e*y
sortsym = x - e*y
n = n_elements(x)

if direction eq 1 then begin
   signx = -1.0 & signy = +1.0
endif else begin
   signx = +1.0 & signy = -1.0
endelse

for i=0,n-1 do begin
    if diff(i) eq 0 then begin
        xi = x(i)
        yi = y(i)   
        sortsymi = sortsym(i)
        j  = n/2
        dj = j
        ;print,'i,x,y,sortsymi=',i,xi,yi,sortsymi
        ;print,'n,j,dj=',n,j,dj
        last = 0
        fail = 0
        while abs(x(j)-signx*xi) gt 1e-6 or abs(y(j)-signy*yi) gt 1e-6 do begin
            if last then begin
                fail = 1
                break
            endif
            
             ;print,'j,dj,sorting=',j,dj,sorting(j)

            if dj eq 1 then last = 1 else dj = (dj+1)/2

            if sortsymi le sorting(j) then j = j - dj else j = j + dj

            if j lt 0 then j = 0
            if j ge n then j = n - 1

        endwhile 

        ;print,'solution j,x,y,sort=',j,x(j),y(j),sorting(j)

        if fail then begin
            print,'error in symdiff: incorrect number of mirror points'
            print,'i,j,dj,last=',i,j,dj,last
            print,'xi,yi=',xi,yi
            if j ge 0 and j lt n then print,'xj,yj=',x(j),y(j)
            diff(i)=1e30
            retall
        endif else begin
            if keyword_set(anti) then diff(i)= a(i)+a(j) else diff(i)= a(i)-a(j)
            diff(j)=-diff(i)
            if report and abs(diff(i)) gt report then print,"i,x,y,diff=",i,xi,yi,diff(i)
        endelse
    endif
endfor
return,diff

end
;===========================================================================
function symmdiffreg,direction,a
;
; Take symmetric difference of "a" with respect to a mirror plane in direction
; "direction"
;
;===========================================================================
on_error,2

if not keyword_set(report) then report = 0

siz=size(a)
dim=siz(0)
nx=siz(1)

diff=a

case dim of
1: for i=0,nx-1 do diff(i)=a(i)-a(nx-1-i)

;if keyword_set(anti) then diff(i)= a(i)+a(nx-1-i) 
;if report and abs(diff(i)) gt report then print,"i,x,y,diff=",i,xi,yi,diff(i)

2: begin
     ny=siz(2)
     case direction of
     1: for i=0,nx-1 do diff(i,*)=a(i,*)-a(nx-1-i,*)
     2: for i=0,ny-1 do diff(*,i)=a(*,i)-a(*,ny-1-i)
     endcase
   end
3: begin
     ny=siz(2)
     nz=siz(3)
     case direction of
     1: for i=0,nx-1 do diff(i,*,*)=a(i,*,*)-a(nx-1-i,*,*)
     2: for i=0,ny-1 do diff(*,i,*)=a(*,i,*)-a(*,ny-1-i,*)
     3: for i=0,nz-1 do diff(*,*,i)=a(*,*,i)-a(*,*,nz-1-i)
     endcase
   end
4: begin
     ny=siz(2)
     nz=siz(3)
     case direction of
     1: for i=0,nx-1 do diff(i,*,*,*)=a(i,*,*,*)-a(nx-1-i,*,*,*)
     2: for i=0,ny-1 do diff(*,i,*,*)=a(*,i,*,*)-a(*,ny-1-i,*,*)
     3: for i=0,nz-1 do diff(*,*,i,*)=a(*,*,i,*)-a(*,*,nz-1-i,*)
     endcase
   end
endcase

return,diff
end

;===========================================================================
function filledge,a

; On the edges use copy of closest cells

siz=size(a)
n1=siz(1)
n2=siz(2)

result=a
result(0,*)   =result(1,*)
result(*,0)   =result(*,1)
result(n1-1,*)=result(n1-2,*)
result(*,n2-1)=result(*,n2-2)

return,result
end

;===========================================================================
pro gengrid,name,x,y,xc,yc,vol2,u,v
;
; From cell center coordinates x,y calculate cell corner coordinates xc,yc,
; cell volumes. Check for array sizes of the optional u,v arguments.
; The name of the calling function is shown for error messages.
;===========================================================================

siz=size(x)
if siz(0) ne 2 then begin
   print,'Function ',name,' is for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

error=''
siz=size(y)
if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd coord'
if keyword_set(u) then begin
   siz=size(u)
   if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='1st func'
endif
if keyword_set(v) then begin
   siz=size(v)
   if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd func'
endif
if error ne '' then begin
  print,'In function ',name,' the first argument does not match the ',error,'.'
  retall
endif

; Coordinates for cell corners
xc=(x(0:n1-2,0:n2-2)+x(0:n1-2,1:n2-1)+x(1:n1-1,0:n2-2)+x(1:n1-1,1:n2-1))/4
yc=(y(0:n1-2,0:n2-2)+y(0:n1-2,1:n2-1)+y(1:n1-1,0:n2-2)+y(1:n1-1,1:n2-1))/4

; Calculate 2*volume=(diagonal_1 X diagonal_2)
vol2=dblarr(n1,n2)+1
vol2(1:n1-2,1:n2-2)= $
 ((xc(1:n1-2,1:n2-2)-xc(0:n1-3,0:n2-3))*(yc(0:n1-3,1:n2-2)-yc(1:n1-2,0:n2-3)) $
 -(yc(1:n1-2,1:n2-2)-yc(0:n1-3,0:n2-3))*(xc(0:n1-3,1:n2-2)-xc(1:n1-2,0:n2-3)))

end

;===========================================================================
function intedge,f,xc
;
; Integrate the neighbouring values of "f" for the four edges described by "xc"
; The size of "f", "xc", and the result are n1*n2, (n1-1)*(n2-1), and n1*n2
; respectively, but only the inner (n1-2)*(n2-2) points are calculated, the
; edge values are 0-s.
;===========================================================================

siz=size(f)
n1=siz(1)
n2=siz(2)

intf=dblarr(n1,n2)
intf(1:n1-2,1:n2-2)=-(xc(1:n1-2,1:n2-2)-xc(0:n1-3,1:n2-2))*f(1:n1-2,2:n2-1) $
                    -(xc(1:n1-2,0:n2-3)-xc(1:n1-2,1:n2-2))*f(2:n1-1,1:n2-2) $
                    -(xc(0:n1-3,0:n2-3)-xc(1:n1-2,0:n2-3))*f(1:n1-2,0:n2-3) $
                    -(xc(0:n1-3,1:n2-2)-xc(0:n1-3,0:n2-3))*f(0:n1-3,1:n2-2)

return,intf

end

;===========================================================================
function intedge_rz,f,rc,zc
;
; Integrate r_edge*f_neighbour*dz for the four cell edges.
; assuming axial symmetry in the ignored direction.
; Only the inner (n1-2)*(n2-2) points are calculated, the edge values are 0-s.
;
;===========================================================================

siz=size(f)
n1=siz(1)
n2=siz(2)

intf=dblarr(n1,n2)
intf(1:n1-2,1:n2-2)= $
    -f(1:n1-2,2:n2-1)*(rc(1:n1-2,1:n2-2)+rc(0:n1-3,1:n2-2)) $
                     *(zc(1:n1-2,1:n2-2)-zc(0:n1-3,1:n2-2)) $
    -f(2:n1-1,1:n2-2)*(rc(1:n1-2,0:n2-3)+rc(1:n1-2,1:n2-2)) $
                     *(zc(1:n1-2,0:n2-3)-zc(1:n1-2,1:n2-2)) $
    -f(1:n1-2,0:n2-3)*(rc(0:n1-3,0:n2-3)+rc(1:n1-2,0:n2-3)) $
                     *(zc(0:n1-3,0:n2-3)-zc(1:n1-2,0:n2-3)) $
    -f(0:n1-3,1:n2-2)*(rc(0:n1-3,1:n2-2)+rc(0:n1-3,0:n2-3)) $
                     *(zc(0:n1-3,1:n2-2)-zc(0:n1-3,0:n2-3))

return,intf

end

;===========================================================================
function grad,idir,f,x,y
;
; Take gradient of "f" in direction "idir" on the "x,y" structured 2D grid.
; Gradient is the contour integral of edge_normal_idir*f_edge_averaged
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Gradient can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(idir) eq 0 or n_elements(f) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function grad'
   retall
endif

gengrid,'grad',x,y,xc,yc,vol2,f

if idir eq 1 then return,filledge( intedge(f,yc)/vol2) $
else              return,filledge(-intedge(f,xc)/vol2)

end

;===========================================================================
function grad_rz,idir,f,r,z
;
; Take gradient of "f" in direction "idir" on the "r,z" structured 2D grid
; assuming axial symmetry in the ignored direction.
; Gradient is the contour integral of edge_normal_idir*f*R_edge_averaged
; divided by R*cell_volume - f/r for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels for idir=2, or equals +f/2R for idir=1.
; Gradient can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(ndir) eq 0 or n_elements(f) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function grad_rz'
   retall
endif

gengrid,'grad_rz',r,z,rc,zc,vol2,f

if idir eq 1 then return,filledge( (intedge_rz(f,rc,zc)/vol2 - f)/2/r ) $
else              return,filledge( -intedge(f,rc^2)/vol2/2/r)

end

;===========================================================================
function div,u,v,x,y
;
; Take divergence of "u,v" vector with respect to "x,y" on a structured 2D grid
; Divergence is the contour integral of edge_normal.(u,v)_edge_averaged
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Divergence can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function div'
   retall
endif

gengrid,'div',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,yc)-intedge(v,xc))/vol2)

end

;===========================================================================
function div_rz,u,v,r,z
;
; Take divergence of "u,v" vector with respect to "r,z" on a structured 2D grid
; assuming axial symmetry in the ignored direction.
; Divergence is the contour integral of edge_normal.(u,v)*R_edge_averaged
; divided by R*cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center is simply u/(2R).
; Divergence can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function div_rz'
   retall
endif

gengrid,'div_rz',r,z,rc,zc,vol2,u,v

return,filledge(((intedge_rz(u,rc,zc)-intedge(v,rc^2))/vol2 + u)/2/r)

end

;===========================================================================
function curl,u,v,x,y
;
; Take curl of "u,v" vector with respect to "x,y" on a structured 2D grid.
; Curl is the contour integral of edge_vector.(u,v)_edge_averaged
; divided by cell_volume for each cells. See also comments for div function.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function curl'
   retall
endif

gengrid,'curl',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,xc)+intedge(v,yc))/vol2)

end

;===========================================================================
function curl_rz,u,v,r,z
;
; Take curl of "u,v" vector with respect to "r,z" on a structured 2D grid
; with axial symmetry in the ignored direction.
; Curl is the contour integral of edge_vector.(u,v)*R_edge_averaged
; divided by R*cell_volume for each cells - v/R.
; See also comments for the div_rz function on edge average and edge cells.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function curl_rz'
   retall
endif

gengrid,'curl',r,z,rc,zc,vol2,u,v

return,filledge(-((intedge_rz(v,rc,zc)+intedge(u,rc^2))/vol2 - v)/2/r)

end

;==========================================
function coarsen,a,boxsize
;
; Produce a coarser array from "a" by averaging out cells in a box.
; The box size can be defined by a scalar (n long interval, n*n square,
; or ,n*n*n cube) or as an
; array of the same dimension as "a" (n1*n2 rectangle or n1*n2*n3 brick)

;on_error,2

if(n_elements(a) eq 0 or n_elements(boxsize) eq 0)then begin
   print,'Calling sequence is: array_co=coarse(array, boxsize)'
   retall
endif

siz=size(a)
ndim=siz(0)

if(ndim eq 0 or ndim gt 3)then begin
   print,'coarse requires a 1,2 or 3D array for the 1st argument'
   retall
endif
nx=siz(1:ndim)

siz=size(box)
if(siz(0) eq 0)then begin
   n=intarr(ndim)+boxsize
endif else if siz(0) eq ndim then begin
   n=boxsize
endif else begin
   print,'boxsize should either be a scalar, or an array '
   print,'of the same dimension as the number of dimensions of the array'
   retall
endelse

case ndim of
   1: begin
      result=dblarr(nx(0)/n(0))
      for ix=0,(nx(0)-1)/n(0) do $
        for i=0,n(0)-1 do $
           result(ix)=result(ix)+a(ix*n(0)+i)
      result=result/n(0)
   end
   2: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1))
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
        for i=0,n(0)-1 do $
        for j=0,n(1)-1 do $
           result(ix,iy)=result(ix,iy)+a(ix*n(0)+i,iy*n(1)+j)
      result=result/n(0)/n(1)
   end
   3: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1),nx(2)/n2)
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
      for iz=0,(nx(2)-1)/n(2) do $
        for i=0,n(0)-1 do $
        for j=0,n(1)-1 do $
        for k=0,n(2)-1 do $
           result(ix,iy,iz)=result(ix,iy,iz)+a(ix*n(0)+i,iy*n(1)+j,iz*n(2)+k)
      result=result/n(0)/n(1)/n(2)
   end
endcase
return,result
end

;===========================================================================
function quadruplet,nx,x0,x1,dx,ny,y0,y1,dy,nz,z0,z1,dz,nw,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=quadruplet(100,0,30,2,100,30,40,1)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(quadruplet(100,0,99,4,100,30,69,2,2,0,1,1))
;===========================================================================

if keyword_set(dx) then begin
   checkdim,1,nx,x0,x1,dx
   all=lindgen(x1+1)
   sub=all(x0:x1)
   ind=sub(where(sub mod dx eq x0 mod dx))
end
if keyword_set(dy) then begin
   checkdim,2,ny,y0,y1,dy
   ixs=ind
   all=lindgen(y1+1)
   sub=all(y0:y1)
   iys=sub(where(sub mod dy eq y0 mod dy))
   ind=(ixs # (0*iys+1)) + ((0*ixs+nx) # iys)
end
if keyword_set(dz) then begin
   checkdim,3,nz,z0,z1,dz
   ixys=ind
   nxy=long(nx)*long(ny)
   all=lindgen(z1+1)
   sub=all(z0:z1)
   izs=sub(where(sub mod dz eq z0 mod dz))
   ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs))
   for iz=0,n_elements(izs)-1 do ind(*,*,iz)=ixys + izs(iz)*nxy
end
if keyword_set(dw) then begin
   checkdim,4,nw,w0,w1,dw
   ixyzs=ind
   nxyz=long(nx)*long(ny)*long(nz)
   all=lindgen(w1+1)
   sub=all(w0:w1)
   iws=sub(where(sub mod dw eq w0 mod dw))
   ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs),n_elements(iws))
   for iw=0,n_elements(iws)-1 do ind(*,*,*,iw)=ixyzs + iws(iw)*nxyz
end

return,ind
end

;===========================================================================
function triplet,x0,x1,dx,y0,y1,dy,z0,z1,dz,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=triplet(0,99,2,0,99,2)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(triplet(0,99,4,0,99,4,0,1,1))
;
; Note: the resulting indices are valid for an array of size
;
;      (x1+1)*(y1+1)*(z1+1)*(w1+1)
;
;===========================================================================

if keyword_set(dw) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz,w1+1,w0,w1,dw)

if keyword_set(dz) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz)

if keyword_set(dy) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy)

if keyword_set(dx) then $
   return,quadruplet(x1+1,x0,x1,dx)

print,'Error in TRIPLET: All strides are 0!'
retall

end

;===========================================================================
pro checkdim,idim,nx,x0,x1,dx

; Check quadruplet for conditions nx>x1>=x0>=0 and dx>0
;===========================================================================

   if nx le 0 then begin
      print,'Size must be positive for dimension',idim
      retall
   endif
   if x1 ge nx then begin
      print,'Maximum index must be less than size for dimension',idim
      retall
   endif
   if x0 lt 0 then begin
      print,'Minimum index must be greater than 0 for dimension',idim
      retall
   endif
   if x0 gt x1 then begin
      print,'Minimum index must be less than maximum index for dimension',idim
      retall
   endif
   if dx le 0 then begin
      print,'Stride must be a positive integer for dimension',idim
      retall
   endif

return
end

;===================================================================
pro plotgrid,x,y,lines=lines,xstyle=xstyle,ystyle=ystyle,polar=polar,$
             xrange=xrange,yrange=yrange,noorigin=noorigin
;===================================================================

on_error,2

if not keyword_set(x) then begin
    print,'Usage: plotgrid, x [,y] [,/lines] [,/polar]',$
      ' [,xstyle=3] [,ystyle=1]',$
      '                   [,xrange=[-10,10]], [yrange=[-10,10]]'
    retall
endif

xx=reform2(x)
sizx=size(xx)

if (n_elements(polar) eq 0) then polar = 0

if not keyword_set(y) then begin
    case sizx(0) of
        3:begin
            if sizx(3) ne 2 then goto, ERROR1
            yy=xx(*,*,1)
            xx=xx(*,*,0)
        end
        2:begin
            if sizx(2) ne 2 then goto, ERROR1
            yy=xx(*,1)
            xx=xx(*,0)
            lines=0
        end
        else: goto, ERROR1
    endcase
endif else begin
    yy=reform2(y)
    sizy=size(yy)
    if sizx(0) ne sizy(0)            then goto, ERROR2
    if max(abs(sizx-sizy)) ne 0      then goto, ERROR2
    if sizx(0) ne 2 and sizx(0) ne 1 then goto, ERROR2
    if sizx(0) eq 1 then lines=0
endelse

if not keyword_set(xrange) then xrange = [0,0]
if not keyword_set(yrange) then yrange = [0,0]

if keyword_set(lines) then begin

    plot, xx, yy, XSTYLE=xstyle, YSTYLE=ystyle, POLAR=polar, $
      XRANGE=xrange, YRANGE=yrange, /NOERASE, /NODATA

    if(keyword_set(noorigin))then begin
        for ix=0,sizx(1)-1 do $
          for iy=0,sizx(2)-2 do $
          if((xx(ix,iy)   ne 0 or yy(ix,iy)   ne 0) and $
             (xx(ix,iy+1) ne 0 or yy(ix,iy+1) ne 0)) then $
          oplot,[xx(ix,iy),xx(ix,iy+1)],[yy(ix,iy),yy(ix,iy+1)],POLAR=polar,$
          psym=0

        for iy=0,sizx(2)-1 do $
          for ix=0,sizx(1)-2 do $
          if((xx(ix,iy)   ne 0 or yy(ix,iy)   ne 0) and $
             (xx(ix+1,iy) ne 0 or yy(ix+1,iy) ne 0)) then $
          oplot,[xx(ix,iy),xx(ix+1,iy)],[yy(ix,iy),yy(ix+1,iy)],POLAR=polar,$
          psym=0

    endif else begin

        for ix=0,sizx(1)-1 do $
          oplot,xx(ix,*),yy(ix,*),POLAR=polar,psym=0
        for iy=0,sizx(2)-1 do $
          oplot,xx(*,iy),yy(*,iy),POLAR=polar,psym=0

    endelse

endif else begin

    if polar then $
      plot, xx, yy, PSYM=3, SYMSIZE=!p.symsize, $
      XRANGE=xrange, YRANGE=yrange, XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE,$
      /POLAR $
    else $
      plot, xx, yy, PSYM=1, SYMSIZE=!p.symsize, $
      XRANGE=xrange, YRANGE=yrange, XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE
endelse

return

ERROR1:
   print,'size(x)=',sizx
   print,'Error: plotgrid,x  requires x(nx,ny,2) array'
   retall

ERROR2:
   print,'size(x)=',sizx,' size(y)=',sizy
   print,'Error: plotgrid,x,y requires x(nx,ny) y(nx,ny) arrays'
   retall


end

;==========================================
pro compare,w0,w1,wnames

; Compare all variables in w0 and w1 by calculating
; relative difference in the 1st norm.
;==========================================

on_error,2

sizew0=size(w0)
sizew1=size(w1)

if sizew0(0) ne sizew1(0) then begin
   print,'w0 and w1 have different dimensions:',sizew0(0),' and ',sizew1(0)
   retall
endif

ndim=sizew0(0)-1

if ndim eq 0 then begin
   ndim=1
   nw=1
endif else $
   nw=sizew0(ndim+1)

if max(abs(sizew0(1:ndim)-sizew1(1:ndim))) gt 0 then begin
   print,'w0 and w1 have different sizes:',sizew0(1:ndim),' /= ',sizew1(1:ndim)
   retall
endif

if keyword_set(wnames) then $
  print,$
  'var 2*max(|A-B|)/max(|A|+|B|) 2*sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)/2' $
else $
  print, $
  'ind 2*max(|A-B|)/max(|A|+|B|) 2*sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)/2'

for iw=0,nw-1 do begin
   case ndim of
   1: begin
      wsum=max(abs(w0(*,iw))+abs(w1(*,iw)))/2
      wdif=max(abs(w0(*,iw)-w1(*,iw)))
      wsum1=total(abs(w0(*,iw))+abs(w1(*,iw)))/2
      wdif1=total(abs(w0(*,iw)-w1(*,iw)))
      end
   2: begin
      wsum=max(abs(w0(*,*,iw))+abs(w1(*,*,iw)))/2
      wdif=max(abs(w0(*,*,iw)-w1(*,*,iw)))
      wsum1=total(abs(w0(*,*,iw))+abs(w1(*,*,iw)))/2
      wdif1=total(abs(w0(*,*,iw)-w1(*,*,iw)))
      end
   3: begin
      wsum=max(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))/2
      wdif=max(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
      wsum1=total(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))/2
      wdif1=total(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
      end
   endcase

   if keyword_set(wnames) then begin
      if wsum eq 0. then print,wnames(iw),' wsum=0' $
      else               print,wnames(iw),wdif/wsum,wdif1/wsum1,wsum
   endif else begin
      if wsum eq 0. then print,iw,' wsum=0' $
      else               print,iw,wdif/wsum,wdif1/wsum1,wsum
   endelse

endfor
end

;=============================================================================
function rel_error, w1, w2, iws

; Calculate relative errors of w1 with respect to w2.
; The 1 to ndim (=1,2, or 3) indexes are the spatial indexes. 
; The last index is assumed to be the variable index.
; If w1 has more elements than w2, it is regarded as the "reference"
; solution, otherwise w2 is regarded as reference solution.
; The difference between w1 and w2 is normalized by the reference solution.
; The function returns the relative errors averaged over all variables.
; If iws is present then the error is calculated for the variable
; indexes listed in iws. 

  if n_elements(w1) le n_elements(w2) then begin
     w     = w1
     wref  = w2
  endif else begin
     w     = w2
     wref  = w1
  endelse

  s     = size(w)
  sref  = size(wref)
  ndim  = s(0) - 1
  nw    = s(ndim+1)

  if nw ne sref(ndim+1) then begin
     print,'ERROR in rel_error: w1 and w2 have different number of variables'
     help,w1,w2
     retall
  endif

  if n_elements(iws) gt 0 then begin
                                ; keep variable indexes iws only in w, wref
     nw = n_elements(iws)
     case ndim of
        1: begin
           w    = w(*,iws)
           wref = wref(*,iws)
        end
        2: begin
           w    = w(*,*,iws)
           wref = wre(*,*,iws)
        end
        3: begin
           w    = w(*,*,*,iws)
           wref = wref(*,*,*,iws)
        end
     endcase
  endif

  ; Coarsen wref if necessary (only in the spatial indexes)
  nx    = s(1:ndim)
  nxref = sref(1:ndim)
  if max(nxref gt nx) then wref = coarsen(wref, [nxref/nx, 1])

  error = 0.0
  for iw = 0, nw-1 do begin
     case ndim of
        1: error = error + (total(abs(w(*,iw)     - wref(*,iw)))) $
                   /        total(abs(wref(*,iw)))
        2: error = error + (total(abs(w(*,*,iw)   - wref(*,*,iw)))) $
                   /        total(abs(wref(*,*,iw)))
        3: error = error + (total(abs(w(*,*,*,iw) - wref(*,*,*,iw)))) $
                   /        total(abs(wref(*,*,*,iw)))
     endcase
  endfor

  return, error/nw

end



;=============================================================================
pro get_log, source, wlog, wlognames, logtime, timeunit, verbose=verbose

; Read the log data from source. If source is an integer, it is 
; interpreted as a unit number. If it is a string, it is taken as the
; filename. Read the content of the file into wlog and the variable 
; names into wlognames. 
; The optional logtime argument is set to the time in hours.
; If verbose is present set show verbose information.
; If versbose is a string, attach it to 'wlog' in the verbose info.

on_error,2

if not keyword_set(source) then begin
   print, $
     'Usage: get_log, source, wlog, wlognames [,logtime, timeunit] [,verbose=verbose]'
   help,source,wlog,wlognames
   retall
endif

itype = size(source,/type)
if itype eq 2 or itype eq 3 then begin
    filesource=0
    unit = source
    file = 'unit '+strtrim(string(unit),2)
    stat = fstat(unit)
    if not stat.open then begin
        print,'get_log error: unit is not open'
        retall
    endif
endif else if itype eq 7 then begin
    filesource=1
    file = source
    unit = 0
    found = 0
    while not found do begin
        unit = unit + 1
        stat = fstat(unit)
        if not stat.open then found = 1
    endwhile
    openr,unit,file
endif else begin
    print,'get_log error: source =',source,$
      ' should be a unit number or a filename.'
    retall
end

if not keyword_set(verbose) then verbose = 0
; If verbose is a string set the index string to it
if size(verbose,/type) eq 7 then index=verbose else index=''

headline=''
readf,unit,headline
wlogname=''
readf,unit,wlogname
i = strpos(wlogname,'#START')
if i ge 0 then strput,wlogname,'      ',i
str2arr,wlogname,wlognames,nwlog

if verbose then begin
  if filesource then print,'logfile',index,'  =',file
  print,'headline',index,' =',strtrim(headline,2)
  for i=0,nwlog-1 do $
    print,FORMAT='("  wlog",A,"(*,",I2,")= ",A)',index,i,wlognames(i)
endif

; Use buffers for efficient reading
buf   = long(10000)
dbuf  = long(10000)
wlog  = dblarr(nwlog,buf)
wlog_ = dblarr(nwlog)
nt    = long(0)
while not eof(unit) do begin
    on_ioerror,close_file
    readf,unit,wlog_
    wlog(*,nt)=wlog_
    nt=nt+1
    if nt ge buf then begin
        buf=buf+dbuf
        wlog=[[wlog],[dblarr(nwlog,buf)]]
    endif
endwhile
close_file: if filesource then close,unit

if verbose then print,'Number of recorded timesteps: nt=',nt

wlog = transpose(wlog(*,0:nt-1))

logtime = log_time(wlog,wlognames,timeunit)
if verbose then print,'Setting logtime',index

end

;=============================================================================

pro plot_log, logfilenames, func, $
              wlog0, wlognames0, wlog1, wlognames1, wlog2, wlognames2,$
              xrange=xrange, yranges=yranges, timeshifts=timeshifts, $
              smooths=smooths, $
              colors=colors, linestyles=linestyles, symbols=symbols, $
              title=title, xtitle=xtitle, ytitles=ytitles, timeunit=timeunit

; Plot variables listed in the space separated func string from the
; files listed in the string array logfilenames.
; Use wlog0...wlognames2 if all the needed arguments are present 
; otherwise read the arrays from the files.
; Use xrange to set the time range, shift times by timeshifts(nlog), 
; set the function ranges with the optional yranges(2,nfunc) array,
; set the colors, the linestyles and the symbols with the 
; colors(nlog), linestyles(nlog) and symbols(nlog) arrays, respectively.
; Smooth data with a boxcar filter of width smooths(nlog). 
; Set the title (default is the list of file names)
; the X title (default is Time [hr]) and 
; the Y titles (defaults are the function names) with the 
; title, xtitle and ytitles(nfunc) arrays, respectively.
; Set the optional variables to zero to get the default behavior.

nlog = n_elements(logfilenames)
str2arr,func,funcs,nfunc

; read in arrays if not present
if   (nlog eq 1 and n_elements(wlognames0) eq 0) $
  or (nlog eq 2 and n_elements(wlognames1) eq 0) $
  or (nlog eq 3 and n_elements(wlognames2) eq 0) then begin

    if nlog ge 1 then get_log,logfilenames(0),wlog0,wlognames0,/verbose
    if nlog ge 2 then get_log,logfilenames(1),wlog1,wlognames1,verbose='1'
    if nlog ge 3 then get_log,logfilenames(2),wlog2,wlognames2,verbose='2'

endif

; Shift times by 0 unless defined
if n_elements(timeshifts) lt nlog then timeshifts = fltarr(nlog)

; Do not smooth data unless defined
if n_elements(smooths) lt nlog then smooths = intarr(nlog)

; Calculate the xrange from the data unless defined
if n_elements(xrange) ne 2 then begin
    DoXrange = 1
    xrange   = [1e30, -1e30]
endif else $
  DoXrange = 0

; Calculate yranges from the data unless defined
if n_elements(yranges) ne 2*nfunc then begin
    DoYrange = 1
    yranges = fltarr(2,nfunc)
    yranges(0,*) =  1e30
    yranges(1,*) = -1e30
endif else $
  DoYrange = 0

; If line styles are not defined use normal lines (0)
if n_elements(linestyles) lt nlog then linestyles = intarr(nlog)

; If symbols are not defined do not use symbols (0)
if n_elements(symbols) lt nlog then symbols = intarr(nlog)

; If colors are not defined, make 3 colors with the last one white/black
if n_elements(colors) lt nlog then colors = intarr(nlog) + 255

; If none of colors, linestyles or symbols are defined, make colors different
if max(linestyles) eq 0 and max(symbols) eq 0 and min(colors) eq 255 then $
  colors = [255,150,250]

; Define default title
if n_elements(title) eq 1 and size(title,/type) eq 7 then $
  title0=title else title0=strjoin(logfilenames,' ')

; Define default xtitle
if n_elements(xtitle) eq 1 and size(xtitle,/type) eq 7 then xtitle0=xtitle $
else if timeunit eq '1' then xtitle0='Time' $
else xtitle0='Time ['+timeunit+']'

; Define default ytitles
if n_elements(ytitles) eq nfunc and size(ytitles,/type) eq 7 then $
  ytitles0=ytitles else ytitles0=funcs

if strpos(!d.name,'X') gt -1 then thick = 1 else thick = 3
if strpos(!d.name,'X') gt -1 then loadct,39 else loadct,40

; Set size of plot windows
ppp   = nfunc
;space = 0.05
;if !y.charsize gt 0 then space = space*!y.charsize
space = 5*float(!d.y_ch_size)/float(!d.y_size)
set_space, ppp, space, sizes, ny = ppp

erase

; The first iteration is used to get the X and Y ranges from the
; data. This can be skipped if both ranges are given explicitly.
if DoXrange or DoYrange then iter0 = 1 else iter0 = 2

for iter = iter0, 2 do begin
    for ilog = 0, nlog-1 do begin
        case ilog of 
            0: begin
                wlog=wlog0
                wlognames = wlognames0
            end
            1: begin
                wlog=wlog1
                wlognames = wlognames1
            end
            2: begin
                wlog=wlog2
                wlognames = wlognames2
            end
        endcase
        
        hour = log_time(wlog,wlognames,timeunit) + timeshifts(ilog)

        for ifunc = 0, nfunc-1 do begin

            field = log_func(wlog, wlognames, funcs(ifunc), error)

            if error then begin
                if iter eq 1 then print,"function ",funcs(ifunc), $
                  " was not found in wlog",strtrim(string(ilog),2)
            endif else begin

                if(smooths(ilog) gt 1)then field = smooth(field,smooths(ilog))

                if iter eq 1 then begin
                    if DoXrange then begin
                        xrange[0]   = min( [ xrange[0], hour ] )
                        xrange[1]   = max( [ xrange[1], hour ] )
                    endif else $
                        field = field( where(hour ge xrange[0] and $
                                             hour le xrange[1]))
                    if DoYrange then begin
                        yranges[0,ifunc] = min( [ yranges[0,ifunc], field ] )
                        yranges[1,ifunc] = max( [ yranges[1,ifunc], field ] )
                    endif
                endif else begin
                    set_position, sizes, 0, ifunc, posm, /rect
                    posm(0) = posm(0) + 0.05
                    if nfunc lt 3 then posm(1) = posm(1) + 0.05/nfunc
                    title1  = ''
                    xtitle1 = ''
                    xtickname1 = strarr(60)+' '
                    xstyle=5
                    ystyle=5
                    if ilog eq 0 then begin
                        xstyle=1
                        ystyle=1
                        if ifunc eq 0       then title1  = title0
                        if ifunc eq nfunc-1 then begin
                            xtitle1    = xtitle0
                            xtickname1 = !x.tickname
                        endif
                    endif
                    plot, hour, field, pos = posm, $
                      xrange = xrange, $
                      yrange = yranges(*,ifunc), $
                      xstyle = xstyle, $
                      ystyle = ystyle, $
                      title  = title1, $
                      xtitle = xtitle1, $
                      xtickname = xtickname1, $
                      ytitle = ytitles0(ifunc), $
                      color = colors(ilog), $
                      psym  = symbols(ilog), $
                      linestyle = linestyles(ilog), $
                      thick = thick, $
                      /noerase

                    if ilog eq nlog-1 then oplot,xrange,[0,0],linestyle=2
                
                endelse
            endelse
        endfor
    endfor
endfor 

if DoXrange then xrange =0
if DoYrange then yranges=0
if strpos(!d.name,'X') lt 0 then loadct,39

end
;=============================================================================
pro rms_logfiles,logfilename,varname,tmin=tmin,tmax=tmax,verbose=verbose

; Print the rms deviation between two logfiles for variables in varname.
; If varname is not present, show rms for all variables.

on_error,2

interpol_logfiles,logfilename,var0,var1,varname,time,tmin=tmin,tmax=tmax,$
  verbose=verbose
str2arr,varname,varnames,nvar
ntime = n_elements(time)

print,'var rms(A-B) rsm(A) rms(B)'
for ivar=0,nvar-1 do $
  print,varnames(ivar),sqrt(total((var0(*,ivar)-var1(*,ivar))^2)/ntime), $
  sqrt(total(var0(*,ivar)^2)/ntime), sqrt(total(var1(*,ivar)^2)/ntime)

end
;============================================================================
pro interpol_logfiles,logfilename,var0,var1,varname,time,tmin=tmin,tmax=tmax,$
                      verbose=verbose

; Interpolate variables between two logfiles for variables in varname.
; If varname is not present, interpolate all variables.

on_error,2

str2arr,logfilename,logfilenames,nfile

get_log, logfilenames(0), wlog0, varnames0, verbose=verbose
get_log, logfilenames(1), wlog1, varnames1, verbose=verbose
if not keyword_set(varname) then varname = varnames0
interpol_log,wlog0,wlog1,var0,var1,varname,varnames0,varnames1,$
  time,tmin=tmin,tmax=tmax

end
;============================================================================
pro interpol_log,wlog0,wlog1,var0,var1,varname,varnames0,varnames1,time,$
             tmin=tmin,tmax=tmax,timeunit=timeunit

; Interpolate the variables listed in varname to the time of wlog0
; between tmin and tmax. 

on_error,2

str2arr,varname,varnames,nvar

if nvar eq 0 then begin
    print,'Usage: interpol_log, wlog0, wlog1, var0, var1, varnames ', $
      '[,varnames0] [,varnames1] [,time] [,tmin=tmin] [,tmax=tmax]'
    retall
endif

if n_elements(varnames0) eq 0 then varnames0 = varnames
if n_elements(varnames1) eq 0 then varnames1 = varnames0

nvar0 = n_elements(varnames0)
nvar1 = n_elements(varnames1)

sizewlog0=size(wlog0)
sizewlog1=size(wlog1)

if sizewlog0(0) ne 2 or sizewlog1(0) ne 2 then begin
   print,'wlog0 and wlog1 must be 2D arrays'
   retall
endif

if sizewlog0(2) ne nvar0 then begin
   print,'Second dimension of wlog0 should be nvar0=',nvar0
   retall
endif

if sizewlog1(2) ne nvar1 then begin
   print,'Second dimension of wlog1 should be nvar1=',nvar1
   retall
endif

time0 = log_time(wlog0,varnames0,timeunit)
time1 = log_time(wlog1,varnames1,timeunit)

if not keyword_set(tmin) then tmin = max([ min(time0), min(time1) ])
if not keyword_set(tmax) then tmax = min([ max(time0), max(time1) ])

index0 = where(time0 ge tmin and time0 le tmax)
time   = time0(index0)
ntime  = n_elements(time)

var0   = fltarr(ntime, nvar)
var1   = fltarr(ntime, nvar)

for ivar = 0, nvar-1 do begin

    name = varnames(ivar)

    field0 = log_func(wlog0, varnames0, name, error0)
    field1 = log_func(wlog1, varnames1, name, error1)

    if error0 then print,'Could not find variable ',name,' in wlog0'
    if error1 then print,'Could not find variable ',name,' in wlog1'

    if not error0 and not error1 then begin

        var0(*,ivar) = field0(index0)
        var1(*,ivar) = interpol(field1,time1,time)

    endif
endfor

end

;==========================================
pro quit
   exit
end
;==========================================

pro set_space, nb, space, sizes, nx = nx, ny = ny

; Determines the size and multiplication factors for plotting perfect circles
; or squares. This routine is used to simply find the parameters, another
; procedure, set_position, is used to actually find the position of the
; circle or square.
; This routine maxamizes the area used by the plots, determinining the best
; positions for the number of plots that the user has selected.
;
; input parameters:
; nb - number of plots on a page
; space - amount of space in between each of the plots in normalized
;          coordinates
;
; output parameters:
; bs - box size (size of the plotting region)
; nbx, nby - number of plots in the x and y directions
; xoff, yoff - x and y offsets for positions
; xf, yf - x and y multiplication factors for making perfect squares
;
; This has been adapted to allow the user to define how many objects
;   are in the x and y direction on Jan 2, 1998

  sizes = {bs:0.0, nbx:0, nby:0, xoff:0.0, yoff:0.0, xf:0.0, yf:0.0, $
           ppp: nb, space:space}

  xsi = float(!d.x_size)
  ysi = float(!d.y_size)

  xs = xsi - 5.0*space*xsi
  ys = ysi - 5.0*space*ysi

  if nb eq 1 then begin

    sizes.nbx = 1
    sizes.nby = 1
    sizes.bs = 1.0 - space

    if xs gt ys then begin

       sizes.yf = 1.0
       sizes.xf = ys/xs

    endif else begin

       sizes.xf = 1.0
       sizes.yf = xs/ys

     endelse

  endif else begin

    if (n_elements(nx) gt 0) then begin
      sizes.nbx = nx(0)
      if n_elements(ny) eq 0 then sizes.nby = nb/nx(0) else sizes.nby = ny(0)
    endif else begin
      if (n_elements(ny) gt 0) then begin
        sizes.nby = ny(0)
        sizes.nbx = nb/ny(0)
      endif else begin
        if xs gt ys then begin
          sizes.nbx = round(sqrt(nb))
          sizes.nby = fix(nb/sizes.nbx)
        endif else begin
          sizes.nby = round(sqrt(nb))
          sizes.nbx = fix(nb/sizes.nby)
        endelse
      endelse
    endelse

    if xs gt ys then begin

      if (sizes.nbx*sizes.nby lt nb) then                               $
        if (sizes.nbx le sizes.nby) then sizes.nbx = sizes.nbx + 1      $
        else sizes.nby = sizes.nby + 1                                  $
      else                                                        	$
	if (sizes.nbx lt sizes.nby) and					$
	   (n_elements(nx) eq 0) and					$
	   (n_elements(ny) eq 0) then begin
	  temp = sizes.nby
	  sizes.nby = sizes.nbx
	  sizes.nbx = temp
	endif

      sizes.yf = 1.0
      sizes.xf = ys/xs
      sizes.bs = ((1.0-space*(sizes.nbx-1))/sizes.nbx )/sizes.xf
      if sizes.nby*sizes.bs+space*(sizes.nby-1) gt 1.0 then 		$
	sizes.bs = (1.0- space*(sizes.nby-1))/sizes.nby

    endif else begin

      if (sizes.nbx*sizes.nby lt nb) then				$
	if (sizes.nby le sizes.nbx) then sizes.nby = sizes.nby + 1	$
	else sizes.nbx = sizes.nbx + 1					$
      else								$
	if (sizes.nby lt sizes.nbx) and					$
	   (n_elements(nx) eq 0) and					$
	   (n_elements(ny) eq 0) then begin
	  temp = sizes.nby
	  sizes.nby = sizes.nbx
	  sizes.nbx = temp
	endif

      sizes.xf = 1.0
      sizes.yf = xs/ys
      sizes.bs = ((1.0 - space*(sizes.nby-1))/sizes.nby)/sizes.yf
      if sizes.nbx*sizes.bs+space*(sizes.nbx-1) gt 1.0 then 		$
	sizes.bs = (1.0 - space*(sizes.nbx-1))/sizes.nbx

    endelse

  endelse

  sizes.xoff = (1.0 - sizes.xf*(sizes.bs*sizes.nbx + space*(sizes.nbx-1)))/2.0
  sizes.yoff = (1.0 - sizes.yf*(sizes.bs*sizes.nby + space*(sizes.nby-1)))/2.0

end

;============================================================================

pro set_position, sizes, xipos, yipos, pos, rect = rect, $
                  xmargin = xmargin, ymargin = ymargin

; used in conjunction with set_space. Determines the position of the current
; plotting region, given the output parameters from set_space.
;
; Input parameters:
; nb, space, bs, nbx, nby, xoff, yoff, xf, yf - Outputs from set_space
; pos_num - the number of the plot, ranges from 0 : bs-1
;
; Output parameters:
;
; pos - the position of the plot, used in the plot command
;
; modified to make rectangles on Jan 2, 1998

  nb = sizes.ppp
  space = sizes.space

  yf2 = sizes.yf
  yf = sizes.yf*(1.0-space)
  xf2 = sizes.xf
  xf = sizes.xf*(1.0-space)

  if keyword_set(rect) then begin

    if keyword_set(xmargin) then xmar = xmargin(0) 			$
    else xmar = space/2.0

    if keyword_set(ymargin) then ymar = ymargin(0) 			$
    else ymar = space/2.0

    xbuffer = 3.0*float(!d.x_ch_size)/float(!d.x_size) * !p.charsize +space/4.0
    xtotal = 1.0 - (space*float(sizes.nbx-1) + xmar + xf2*space/2.0) - xbuffer
    xbs = xtotal/(float(sizes.nbx)*xf)

    xoff = xmar - xf2*space/2.0 + xbuffer - space/4.0

    ybuffer = 3.0*float(!d.y_ch_size)/float(!d.y_size) * !p.charsize
    ytotal = 1.0 - (space*float(sizes.nby-1) + ymar + yf2*space/2.0) - ybuffer
    ybs = ytotal/(float(sizes.nby)*yf)

    yoff = space/4.0

  endif else begin

    xbs  = sizes.bs
    xoff = sizes.xoff
    ybs  = sizes.bs
    yoff = sizes.yoff

  endelse

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0 + xbs*xf

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space + xbs*xf

  ypos0 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf + ybs*yf) - yoff
  ypos1 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf) - yoff

  pos= [xpos0,ypos0,xpos1,ypos1]

end

;============================================================================
pro plot_color_bar, pos, maxmin

; plot color bar based on the current color table

xrange=!x.range & yrange=!y.range & !x.range=0 & !y.range=0

maxi = max(maxmin)
mini = min(maxmin)

array = findgen(10,256)
for i=0,9 do array(i,*) = findgen(256)/(256-1)*(maxi-mini) + mini

levels=(findgen(60)-1)/(58-1)*(maxi-mini)+mini

contour, array, /noerase, /cell_fill, xstyle = 5, ystyle = 5, $
  levels = levels, pos=pos, title=' '

plot, maxmin, /noerase, pos = pos, xstyle=1, ystyle=1, /nodata,$
  xtickname = [' ',' '], xticks = 1, xminor=1  , $
  ytickname = strarr(60) + ' ', yticklen = 0.25, $
  title=' ', xtitle=' ',ytitle=' '
axis, 1, ystyle=1, /nodata, yax=1, charsize=0.9*(!p.charsize > 1.), $
  ytitle=' ', ytickname = strarr(60)

!x.range=xrange & !y.range=yrange

end

;============================================================================
pro makect, color

; Create color table corresponding to color='mid','blue','red','rwb','bwr'

common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

; Get number of colors
n=!d.table_size
if n lt 10 or n gt 256 then n=256

r = fltarr(n)
g = fltarr(n)
b = fltarr(n)

if not keyword_set(color) then begin

    print,'red   - white to red'
    print,'blue  - white to blue'
    print,'rwb   - red white blue'
    print,'bwr   - blue white red'
    print,'mid   - blue green white yellow red'

    color = ''
    read,'Enter color table from list above : ', color

endif

color = strlowcase(color)

; Set read, green, blue to values normalized to the 0.0 -- 1.0 range.

case color of
    'red' : begin
        r(*) = 1.
        g(*) = 1. - findgen(n)/(n-1)
        b(*) = 1. - findgen(n)/(n-1)
    end

    'blue' : begin
        r(*) = 1. - findgen(n)/(n-1)
        b(*) = 1.
        g(*) = 1. - findgen(n)/(n-1)
    end

    'rwb' : begin
        half=n/2
        r(0:half-1) = 1.
        g(0:half-1) = findgen(half)/(half-1)
        b(0:half-1) = findgen(half)/(half-1)

        r(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        b(half:n-1) = 1.
    end

    'bwr' : begin
        half=n/2
        b(0:half-1) = 1.
        g(0:half-1) = findgen(half)/(half-1)
        r(0:half-1) = findgen(half)/(half-1)
        
        b(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        r(half:n-1) = 1.
    end

    'mid' : begin
        r(0:n/3-1)     = 0.0
        r(n/3:n/2-1)   = findgen(n/2-n/3)/(n/2-n/3-1)
        r(n/2:n-1)     = 1.0
        
        b(0:n/2-1)      = 1.
        b(n/2:2*n/3-1)  = 1. - findgen(2*n/3-n/2)/(2*n/3-n/2-1)
        b(2*n/3-1:n-1)  = 0.
        
        g(0:n/3-1)      = findgen(n/3)/(n/3-1)
        g(n/3:2*n/3-1)  = 1.
        g(2*n/3:n-1)    = 1. - findgen(n-2*n/3)/(n-2*n/3-1)
        
    end

    else : begin
        print, "Unknown value for color=",color
        r(*) = findgen(n)
        g(*) = findgen(n)
        b(*) = findgen(n)
    end

endcase

r(0) = 0.0
g(0) = 0.0
b(0) = 0.0

r(n-1) = 1.0
g(n-1) = 1.0
b(n-1) = 1.0

r=255*r
g=255*g
b=255*b

r_orig = r
g_orig = g
b_orig = b
r_curr = r_orig
g_curr = g_orig
b_curr = b_orig
tvlct,r,g,b

end
