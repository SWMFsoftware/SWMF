;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
; this program is not as well written as it could be

;********************************
function yticks, axis,index,value
power=alog10(value)
return, string(power, format="('10!a',i2,'!n')")
end

;************
pro psbeg,q1,filename,xsz,ysz

q1 = strupcase(q1)

if ((q1 ne 'Y') and (q1 ne 'N')) then begin
 
 qa1:
 print, 'Make PostScript file ? (Y/N)'
 q1=' '
 read,q1 

 q1 = strupcase(q1)

 if ((q1 ne 'Y') and (q1 ne 'N')) then begin
   goto, qa1
 endif

endif

if (q1 eq 'Y') then begin
  set_plot, 'ps'
  device, file=filename+'.ps'
  wide=xsz*.01
  tall=ysz*.01
  xshift=(8.5-wide)/2.0
  yshift=(11.25-tall)/2.0
  device,$
    color=1,bit=8,/port,xsi=wide,ysi=tall,xoff=xshift, $
  yoff=yshift, /inch, /times
endif else begin
  set_plot,'X'
;  device,retain=2
 window,0,title='Plasmasphere',xsize=xsz,ysize=ysz
endelse

end

;*********
pro psend

device, /close
set_plot,'x'

end
;******
pro clr,col,rev

loadct,col,/silent

tvlct,r,g,b,/get

if rev eq -1 then begin
   r = reverse(r)
   g = reverse(g)
   b = reverse(b)
endif

ia=!d.table_size-1
;print,size(ia)
print,ia

ib=1

r[0:ib] = ia  ;255
g[0:ib] = ia  ;255
b[0:ib] = ia  ;255
r[ib+1] = 0
g[ib+1] = 0
b[ib+1] = 0

tvlct,r,g,b

!p.background = 0
!p.color = ib+1
erase,ia

end

;;;;;;;;;;;;;;;;
;;; pro open ;;;
;;;;;;;;;;;;;;;;

   pro open,filename,data

   on_ioerror,err1

   openr,1,filename

;vphicells,vthetacells (spherical coordinates) in degrees
;of ionospheric foot points,
;theta is zero at the pole and phi is zero at 24 MLT positive towards dawn.
;mgridx, mgridy are GSM equatorial crossing points in Re.

   nthetacells = 0
   nphicells = 0
   readf,1,nthetacells,nphicells
   vthetacells = fltarr(nthetacells)
   vphicells = fltarr(nphicells)
   readf,1,vthetacells
   readf,1,vphicells
   mgridden = fltarr(nthetacells,nphicells)
   mgridx = fltarr(nthetacells,nphicells)
   mgridy = fltarr(nthetacells,nphicells)
   mgridoc = fltarr(nthetacells,nphicells)
   readf,1,mgridden
   readf,1,mgridx
   readf,1,mgridy
   readf,1,mgridoc

   close,1

   print,'Read file '+filename

   datastruct = {nthetacells : 0.0,                             $
                   nphicells : 0.0,                             $
                 vthetacells : fltarr(nthetacells),             $
                   vphicells : fltarr(nphicells),               $
                    mgridden : fltarr(nthetacells,nphicells),   $
                      mgridx : fltarr(nthetacells,nphicells),   $
                      mgridy : fltarr(nthetacells,nphicells),   $
                     mgridoc : fltarr(nthetacells,nphicells)   }

   data = ptr_new(datastruct)

   (*data).nthetacells = nthetacells
   (*data).nphicells = nphicells
   (*data).vthetacells = vthetacells
   (*data).vphicells = vphicells
   (*data).mgridden = mgridden
; the axis are flipped around here so that when plotted +x is up, and +y is to the left
   (*data).mgridx = -mgridy  
   (*data).mgridy = mgridx
   (*data).mgridoc = mgridoc

   return

   err1:
   close,1
   print,'Error reading file '+filename
   return

   end

;;;;;;;;;;;;;;;;
;;; makepoly ;;;
;;;;;;;;;;;;;;;;

   pro makepoly, q1, data, dzrange,xp,yp,xs,ys,xsz,ysz

   rad = !pi / 180.0      

   mgridden = (*data).mgridden
   mgridx = (*data).mgridx
   mgridy = (*data).mgridy
   mgridoc = (*data).mgridoc
   nthetacells = (*data).nthetacells
   nphicells = (*data).nphicells
   vthetacells = (*data).vthetacells
   vphicells = (*data).vphicells

   r=sqrt(mgridx*mgridx + mgridy*mgridy)
   modden=mgridden*(r^4.0)/max(mgridden)

   ;xmx=18.0
   ;xmn=-18.0
   ;ymx=13.0
   ;ymn=-15.0

;   xmx=15.1
;   xmn=-15.1
;   ymx=15.1
;   ymn=-15.1

   xmx=10.1
   xmn=-10.1
   ymx=10.1
   ymn=-10.1

   szx=xs
   szy=ys*((ymx-ymn)/(xmx-xmn))

   plot,[xmx,ymx,-xmn,-ymn],[0,1.57,3.1415,4.71],$
     /nodata,$
     position=[xp,yp,xp+szx,yp+szy],$
     /nor,/noerase,xsty=5,ysty=5,/polar

   print, 'Negative values: ',where(mgridden lt 0.0);,mgridden(where(mgridden lt 0.0))

   logden=mgridden
   logden(where(mgridoc eq 1.0)) = alog10(1e-6*mgridden(where(mgridoc eq 1.0)))

   print,'Triangulating'
   triangulate,mgridx,mgridy,triangles
   print,'Done Triangulating'
   delx=0.02
   dely=0.02
   trigridden = trigrid(mgridx,mgridy,logden,triangles,[delx,dely],[xmn,ymn,xmx,ymx])

   xgridno = ((xmx - xmn)/delx) + 1
   ygridno = ((ymx - ymn)/dely) + 1

   trigridx = findgen(xgridno)*delx + xmn
   trigridy = findgen(ygridno)*dely + ymn

   xgrid=trigridx # REPLICATE(1.,ygridno)
   ygrid=TRANSPOSE( trigridy # REPLICATE(1.,xgridno) )

   r = sqrt(xgrid*xgrid + ygrid*ygrid)

   trigridden(where(r lt 1.0)) = 0.0

   if (q1 eq 'N') then begin
      congridden = congrid(trigridden,szx*xsz,szy*ysz)
	print, 'zrange=',min(congridden),max(congridden)
      denscl = 1b + bytscl(congridden,top=!d.n_colors-3,min=dzrange[0],max=dzrange[1])
;      denscl = 1b + bytscl(congridden,top=!d.table_size-3,min=dzrange[0],max=dzrange[1])
	print,'denscl range:',min(denscl),max(denscl)
      tv,reverse(transpose(denscl)),xp,yp,xsize=szx,ysize=szy,/nor
   endif else begin
      denscl = 1b + bytscl(trigridden,top=!d.n_colors-3,min=dzrange[0],max=dzrange[1])
      tv,reverse(transpose(denscl)),xp,yp,xsize=szx,ysize=szy,/nor
   endelse

   axis,xax=0,0,0,xsty=1,xtickname=['10','5',' ','5','10']
   axis,yax=0,0,0,ysty=1,ytickname=['10','5',' ','5','10']

   end

;*********************
pro colorbar,xp,yp,xs,ys,t,s,xsz,ysz

scale=((findgen(7)/6.0)*(s-t))+t

sden=fltarr(2,7)
sden(0,*)=scale
sden(1,*)=scale
sden1=tri_surf(sden,nx=300,ny=300)

szx=xs
szy=ys

;picto=1b+bytscl(sden1,min=t,max=s,top=!d.n_colors-3)
picto=1b+bytscl(sden1,min=t,max=s,top=!d.table_size-3)
print, 'Bar scale:',min(picto),max(picto)
;print,transpose(sden1(0,*))
print,transpose(picto(0,*))
tv,congrid(picto,szx*xsz,szy*ysz),xp,yp,xsize=szx,ysize=szy,/nor
plot,[0,1],[10.^(t),10.^(s)],pos=[xp,yp,xp+szx,yp+szy],$
    /nor,xsty=5,ysty=5,$
    /noerase,yty=1,/nodata
axis,yaxis=1,ysty=1,yty=1,ticklen=-.5,ytickform='yticks',charsize=1.25,charthick=1.9
axis,yaxis=0,ytitle='Density ( cm!u-3!n )',yticks=1,yticklen=0,$
    ytickn=[' ',' '],charsize=1.75,charthick=1.9

end

;************
;
;   next1.pro
;
;************

pro eqden,filename,q1

;filename = 'plasmasphere025_1

If n_params(0) lt 2 then q1=''

xsz=650.
plotsize=xsz-200.
ysz=plotsize+80.
xplace=20./xsz
yplace=20./ysz
xbarplace=(100.+plotsize)/xsz
xplot=plotsize/xsz
yplot=plotsize/ysz
xbar=30./xsz
ybar=1.*yplot
yhour=(ysz-45.)/ysz

psbeg,q1,filename,xsz,ysz
clr,39,1

data = ptr_new()
open,filename+'.dat',data

dzrange=[0.0,4.0]

makepoly,q1,data,dzrange, xplace, yplace, xplot, yplot, xsz, ysz

colorbar,xbarplace,yplace,xbar,ybar,dzrange[0],dzrange[1],xsz,ysz

;printhour,filename,yhour

if (q1 eq 'Y') then psend

print,'finished'
end
