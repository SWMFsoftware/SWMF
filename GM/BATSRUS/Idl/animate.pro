;^CFG COPYRIGHT VAC_UM
;===========================================================================
;    Written by G. Toth for the Versatile Advection Code.
;
;    Read pictures from one or more ascii or binary ini or out file and 
;    plot or animate one or more functions of w using different plotting 
;    routines. The functions are defined in the "Idl/funcdef.pro" file.
;
;    For generalized coordinates the variables can be interpolated from the 
;    irregular grid onto a regular one, or to polar coordinates.
;
;    A subset can be cut from the grid by using the "cut" index array, e.g.:
;    cut=grid(10:30,*), where "grid" contains the full index array.
;    for the regular grid. The grid array is only defined after animate ran.
;
;    Usage:
;
; .r animate
;
;===========================================================================

   print,'======= CURRENT ANIMATION PARAMETERS ================'
   print,'firstpict=',firstpict,', dpict=',dpict,', npictmax=',npictmax, $
     FORMAT='(a,'+string(n_elements(firstpict))+'i4,a,' $
     +string(n_elements(dpict))+'i4,a,i4)'
   print,'savemovie (n/ps/png/tiff/bmp/jpeg)=',savemovie
   print,'ax,az=',ax,',',az,', contourlevel=',contourlevel,$
     ', velvector=',velvector,', velspeed (0..5)=',velspeed,$
     FORMAT='(a,i4,a,i3,a,i3,a,i4,a,i2)'
   if keyword_set(multiplot) then begin
        siz=size(multiplot)
        ; scalar multiplot value is converted to a row (+) or a column (-)
        if siz(0) eq 0 then begin
            if multiplot gt 0 then multiplot=[multiplot,1,1] $
            else                   multiplot=[1,-multiplot,1]
        endif
        print,'multiplot= ',multiplot,', axistype (coord/cells)=',axistype,$
              ', fixaspect= ',fixaspect,$
              FORMAT='(a,"[",i2,",",i2,",",i2,"]",a,a,a,i1)'
   endif else $
        print,'multiplot= 0 (default), axistype (coord/cells)=',axistype,$
              ', fixaspect= ',fixaspect,$
              FORMAT='(a,a,a,i1)'
   print,'bottomline=',bottomline,', headerline=',headerline,$
              FORMAT='(a,i1,a,i1)'

   if keyword_set(cut) then help,cut
   if keyword_set(wsubtract) then help,wsubtract
   if keyword_set(velpos) then help,velpos
   velpos0 = velpos

   print,'======= FILE DESCRIPTION ============================'
   nfile=0
   if filename eq '' and logfilename ne '' then begin
      filename=logfilename
      while strpos(filename,'.log') ge 0 $
         do strput,filename,'.out',strpos(filename,'.log')
      askstr,'filename(s)   ',filename,1
   endif else $
      askstr,'filename(s)   ',filename,doask

   if stregex(filename, '[?*[]', /boolean) then begin
       spawn,'/bin/ls '+filename, filenames
       nfile = n_elements(filenames)
   endif else $
     str2arr,filename,filenames,nfile
   gettype,filenames,filetypes,npictinfiles
   print,'filetype(s)   =',filetypes
   print,'npictinfile(s)=',npictinfiles

   ; Extend firstpict and dpict into arrays of size nfile
   arr2arr,firstpict,nfile
   arr2arr,dpict,nfile

   ;====== OPEN FILE(S) AND READ AND PRINT HEADER(S)

   anygencoord=0
   for ifile=0,nfile-1 do begin
      openfile,10,filenames(ifile),filetypes(ifile)
      gethead,10,filetypes(ifile), $
         headline,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables
      anygencoord=anygencoord or gencoord
      print,         'headline                  =',strtrim(headline,2)
      print,FORMAT='("variables                 =",100(a," "),$)',variables
      print,FORMAT='(" (ndim=",i2,", nw=",i2,")")',ndim,nw
   endfor

   print,'======= PLOTTING PARAMETERS ========================='
   readplotpar,ndim,cut,cut0,plotdim,nfunc,func,funcs,funcs1,funcs2,$
      nplot,plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask

   readtransform,ndim,nx,anygencoord,transform,nxreg,xreglimits,wregpad,$
     nvector,vectors,grid,doask

   print,'======= DETERMINE PLOTTING RANGES ==================='

   readlimits,nfunc,funcs,autoranges,noautorange,fmax,fmin,doask

   if noautorange then begin
      npict = min( (npictinfiles-firstpict)/dpict + 1 )
      if npict gt npictmax then npict=npictmax
      if npict lt 0 then npict=0
   endif else begin
      npict=0
      for ifile=0,nfile-1 do $
         openfile,ifile+10,filenames(ifile),filetypes(ifile)
      error=0
      while npict lt npictmax and not error do begin

         for ifile=0,nfile-1 do begin

            if npict eq 0 then nextpict=firstpict(ifile) $
            else               nextpict=dpict(ifile)

            get_pict,ifile+10,filetypes(ifile),nextpict,x,w,$
                headline, it, time, gencoord, ndim, neqpar, nw, nx,$
                eqpar, variables, rBody, err

            if keyword_set(wsubtract) then w=w-wsubtract
            wnames=variables(ndim:ndim+nw-1)
            error=err or error

            if not error then begin
                do_transform,transform,ifile,gencoord,variables,nw,x,w, $
                  xreg,wreg,nxreg,xreglimits,x_old,nxreg_old,xreglimits_old,$
                  wregpad,triangles,symmtri,nvector,vectors,usereg

               first= npict eq 0 and ifile eq 0
               getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,fmax,fmin,$
                 doask,x,w,xreg,wreg,usereg,time,eqpar,variables,$
                 cut0,rcut

               if ifile eq nfile-1 then begin
                  if npict eq 0 then print,FORMAT='("ipict:    ",$)'
                  npict=npict+1
                  print,FORMAT='(i4,$)',npict
               endif
            endif
         endfor
      endwhile
      print
      for ifunc=0,nfunc-1 do $
      print,'Min and max value for ',funcs(ifunc),':',fmin(ifunc),fmax(ifunc)

   endelse
   print,'npict=',npict
   if npict eq 0 then begin
      print,'There are no frames to animate! Check the following settings:'
      print,'   npictinfiles=',npictinfiles
      print,'   firstpict   =',firstpict
      print,'   dpict       =',dpict
      print,'   npictmax    =',npictmax
      if min(npictinfiles - firstpict) lt 0 then $
        print,'   firstpict is larger than npictinfiles for some files!' 
      retall
   endif

   ;===== DO ANIMATION IN MULTIX * MULTIY MULTIPLE WINDOWS

   if keyword_set(multiplot) then begin
      multix=multiplot(0)
      multiy=multiplot(1)
      multidir=multiplot(2)
      npict1=(multix*multiy)/(nplot*nfile)
      if npict1 eq 0 then npict1=1
   endif else if nfile eq 1 then begin
      multix=long(sqrt(nplot-1)+1)
      multiy=long((nplot-1)/multix+1)
      multidir=0
      npict1=1
   endif else begin
      multix=nfile
      multiy=nplot
      multidir=1
      npict1=1
   endelse


   if savemovie ne 'n' then spawn,'/bin/mkdir -p Movie'
   if savemovie eq 'ps' then set_plot,'PS',/INTERPOLATE

   doanimate= npict gt npict1 and !d.name eq 'X'
   if !d.name eq 'X' then begin
       if !d.window lt 0 then window
       wshow
   endif
   if doanimate then xinteranimate,set=[!d.x_size,!d.y_size,(npict-1)/npict1+1]

   ipict=0
   ipict1=0
   iplot=0
   for ifile=0,nfile-1 do openfile,ifile+10,filenames(ifile),filetypes(ifile)
   error=0
   while ipict lt npict and not error do begin
      if ipict1 eq 0 then begin
         if not keyword_set(noerase) then erase
         !p.multi=[0,multix,multiy,0,multidir]
         if savemovie eq 'ps' then $
            device,filename='Movie/'+string(FORMAT='(i4.4)',iplot+1)+'.ps',$
		   XSIZE=24,YSIZE=18,/LANDSCAPE,/COLOR,BITS=8
      endif

      if ipict eq 0 then print,FORMAT='("ipict:    ",$)'
      print,FORMAT='(i4,$)',ipict+1

      for ifile=0,nfile-1 do begin

         if ipict eq 0 then nextpict=firstpict(ifile) $
         else               nextpict=dpict(ifile)

         if npict gt 1 or nfile gt 1 or noautorange then begin
            get_pict, ifile+10, filetypes(ifile), nextpict, x, w, $
               headline, it, time, gencoord, ndim, neqpar, nw, nx,$
               eqpar, variables, rBody, err

            error=error or err
         endif

         if not error then begin

            if keyword_set(wsubtract) then w=w-wsubtract

            wnames=variables(ndim:ndim+nw-1)

            do_transform,transform,ifile,gencoord,variables,nw,x,w, $
              xreg,wreg,nxreg,xreglimits,x_old,nxreg_old,xreglimits_old,$
              wregpad,triangles,symmtri,nvector,vectors,usereg

	    linestyle=0
            if multix*multiy lt nplot*nfile then $
              linestyle=fix(nplot*ifile)/(multix*multiy)

            if(keyword_set(timetitle))then begin
                t = time
                if(keyword_set(timetitleunit))then  t = t / timetitleunit
                if(keyword_set(timetitlestart))then t = t - timetitlestart
                plottitles(*) = string(format=timetitle,t)
            endif

            plot_func,x,w,xreg,wreg,usereg,ndim,time,eqpar,rBody,$
              variables,wnames,axistype,plotmodes,plottitles,$
              ax,az,contourlevel,linestyle,$
	      velvector,velspeed,velseed,velpos,velx,vely,veltri,$
              cut,cut0,rcut,plotdim,$
              nfunc,multix,multiy,fixaspect,plotix,plotiy,$
              funcs,funcs1,funcs2,fmin,fmax,f

            if npict1 le 1 then begin
               putbottom,multix,multiy,ifile,0,bottomline,nx,it,time
               putheader,nfile,1,ifile,0,headerline,headline,nx
            endif else begin
               putbottom,multix,multiy,plotix,plotiy,bottomline,nx,it,time
               if ipict1 eq 0 then $
                  putheader,nfile,1,ifile,0,headerline,headline,nx
            endelse
         endif

      endfor

      if ipict1 eq npict1-1 or ipict eq npict-1 then begin
         if doanimate then $
            xinteranimate,frame=iplot,window=!d.window
         if savemovie eq 'ps' then begin
            print,FORMAT='(" (Movie/",i4.4,".ps)",$)',iplot+1
            device,/close
         endif else if savemovie ne 'n' and !d.name eq 'X' then begin
            imagefile=string(FORMAT='("Movie/",i4.4,".",a)',iplot+1,savemovie)
            print,FORMAT='("(",a,")",$)',imagefile
            common colors,  r_curr, g_curr, b_curr
	    write_image,imagefile,savemovie, $
              tvrd( order=(savemovie eq 'tiff'), true=1),r_curr, g_curr, b_curr
         endif
      endif

      ipict1=ipict1+1
      if ipict1 ge npict1 then begin
         ipict1=0
         iplot=iplot+1
      endif
      ipict=ipict+1

   endwhile

   for ifile=0,nfile-1 do close,ifile+10
   print
   !p.multi=0
   !p.title=''
   !x.title=''
   !y.title=''
   !z.title=''
   if savemovie eq 'ps' then set_plot,'X'
   ; Restore velpos array
   velpos=velpos0 & velpos0=0
   if doanimate then xinteranimate,5,/keep_pixmaps
end
