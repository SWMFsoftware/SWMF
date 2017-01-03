;  Copyright (C) 2002 Regents of the University of Michigan, 
;  portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;===========================================================================
;    Written by G. Toth for the Versatile Advection Code
;
;    Animate slices of 3D data read by getpict or animate.
;    The functions are defined in the "Idl/funcdef.pro" file.
;    "cut" should be defined in 2D, possibly using the grid2d index array.
;
;    Usage:
;
; .r slice
;
;===========================================================================

   if not keyword_set(nfile) then begin
      print,'No file has been read yet, run getpict or animate!'
      return
   endif

   if nfile gt 1 then begin
      print,'More than one files were read...'
      print,'Probably w is from file ',filenames(nfile-1)
      nfile=1
   endif

   siz=size(x)

   if siz(0) ne 4 then begin
      print,'The x array is not 3 dimensional'
      return
   endif
   n1=siz(1) & n2=siz(2) & n3=siz(3)

   print,'======= CURRENT SLICING PARAMETERS ================'
   print,'firstslice=',firstslice,', dslice=',dslice,$
        ', nslicemax=',nslicemax,', savemovie (y/n)=',savemovie,$
        FORMAT='(a,i4,a,i4,a,i4,a,a)'
   print,'ax,az=',ax,',',az,', contourlevel=',contourlevel,$
         ', velvector=',velvector,', velspeed (0..5)=',velspeed,$
        FORMAT='(a,i4,a,i3,a,i3,a,i4,a,i2)'
   if keyword_set(multiplot) then begin
        siz=size(multiplot)
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
   if keyword_set(velpos) then help,velpos
   velpos0=velpos

   print
   help,x
   asknum,'slicedir (1, 2, or 3)',slicedir,doask
   asknum,'firstslice',firstslice,doask
   asknum,'dslice',dslice,doask

   siz=size(x)
   if dslice gt 0 then nslice=(siz(slicedir)-firstslice)/dslice+1 $
   else                nslice=-firstslice/dslice+1
   if nslice gt nslicemax then nslice=nslicemax
   if nslice lt 1 then nslice=0
   print,'Number of slices:',nslice
   if nslice lt 1 then begin
       print,'There are no slices to be shown!'
       print,'   Check firstslice=',firstslice,' and dslice=',dslice
       if siz(slicedir) lt firstslice then $
       print,'   The value of firstslice is larger than the'
       print,'   grid size in slice direction (',slicedir,') =',siz(slicedir)
       retall
   endif

   ; store 3D data
   var3d = variables
   x3d   = x
   w3d   = w

   case slicedir of
        1:begin
           x=dblarr(n2,n3,2)
           w=dblarr(n2,n3,nw)
           variables=var3d(1:*)
           grid2d=lindgen(n2,n3)
        end
        2:begin
           x=dblarr(n1,n3,2)
           w=dblarr(n1,n3,nw)
           variables=[var3d(0),var3d(2:*)]
           grid2d=lindgen(n1,n3)
        end
        3:begin
           x=dblarr(n1,n2,2)
           w=dblarr(n1,n2,nw)
           variables=[var3d(0:1),var3d(3:*)]
           grid2d=lindgen(n1,n2)
        end
   endcase

   help,grid2d

   print,'======= PLOTTING PARAMETERS ========================='
   read_plot_param

   usereg=0

   if keyword_set(multiplot) then begin
      multix=multiplot(0)
      multiy=multiplot(1)
      nslice1=(multix*multiy)/(nplot*nfile)
      !p.multi=[0,multix,multiy,0,multiplot(2)]
   endif else begin
      multix=long(sqrt(nplot-1)+1)
      multiy=long((nplot-1)/multix+1)
      nslice1=1
      !p.multi=[0,multix,multiy,0,0]
   endelse

   print,'======= DETERMINE PLOTTING RANGES ==================='
   read_limits

   if not noautorange then $
   for islice=1,nslice do begin
      ix=dslice*(islice-1)+firstslice-1
      case slicedir of
        1:begin
           x(*,*,*)=x3d(ix,*,*,1:2)
           w(*,*,*)=w3d(ix,*,*,*)
        end
        2:begin
           x(*,*,0)=x3d(*,ix,*,0)
           x(*,*,1)=x3d(*,ix,*,2)
           w(*,*,*)=w3d(*,ix,*,*)
        end
        3:begin
           x(*,*,*)=x3d(*,*,ix,0:1)
           w(*,*,*)=w3d(*,*,ix,*)
        end
     endcase
    
      first= islice eq 1
      get_limits,first

   endfor

   print
   for ifunc=0,nfunc-1 do $
        print,'Min and max value for ',funcs(ifunc),':',fmin(ifunc),fmax(ifunc)


   ;==== ANIMATE THE SLICES
   doanimate= nslice gt nslice1 and !d.name eq 'X'
   if !d.name eq 'X' then begin
       if !d.window lt 0 then window
       wshow
   endif
   if doanimate then xinteranimate,set=[!d.x_size,!d.y_size,nslice]

   islice1=0 ; slice index in a multiplot frame
   iplot=0   ; plot index for animation
   for islice=1,nslice do begin
      ix=dslice*(islice-1)+firstslice-1
      dirname=variables(slicedir-1)
      case slicedir of
        1:begin
           x(*,*,*)=x3d(ix,*,*,1:2)
           w(*,*,*)=w3d(ix,*,*,*)
           height=x3d(ix,0,0,0)
           info1='i'+dirname+'='+string(ix,format='(i4)')
        end
        2:begin
           x(*,*,0)=x3d(*,ix,*,0)
           x(*,*,1)=x3d(*,ix,*,2)
           w(*,*,*)=w3d(*,ix,*,*)
           height=x3d(0,ix,0,1)
           info1='i'+dirname+'='+string(ix,format='(i4)')
        end
        3:begin
           x(*,*,*)=x3d(*,*,ix,0:1)
           w(*,*,*)=w3d(*,*,ix,*)
           height=x3d(0,0,ix,2)
           info1='i'+dirname+'='+string(ix,format='(i4)')
        end
      endcase
      info2=dirname+'='+string(height)
      case bottomline of
      0: info=''
      1: info=info2
      2: info=info1+' '+info2
      3: info='time='+string(time,format='(g12.5)')+' '+info1+' '+info2
      endcase

      if not keyword_set(noerase) and islice1 eq 0 then erase

      if velrandom then velpos=0

      rBody3d = rBody
      if abs(height) ge rBody3d then rBody=0.0 $
      else                           rBody=sqrt(rBody3d^2 - height^2)

      plot_func

      xyouts,5+(plotix*!d.x_size)/multix, $
             (1-dyslicelabel)*!d.y_size $
             +  dyslicelabel*(plotiy*!d.y_size)/multiy,/DEV,info
      putheader,1,1,0,0,headerline,headline,nx

      if doanimate then xinteranimate,frame=iplot,window=!d.window
   
      islice1=islice1+1
      if islice1 ge nslice1 then begin
         islice1=0
         iplot=iplot+1
      endif
   endfor

   ; restore 3d state
   x         = x3d
   w         = w3d
   variables = var3d
   rBody     = rBody3d
   
   print
   !p.multi=0
   !p.title=''
   !x.title=''
   !y.title=''
   !z.title=''
   ; Restore velpos array
   velpos=velpos0 & velpos0=0

   if doanimate then xinteranimate,5,/keep_pixmaps
end
