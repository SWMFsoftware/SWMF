; READ_GEOMETRY.PRO

; bdir

  bdirsx = fltarr(nz,nf,nl)
  bdirsy = fltarr(nz,nf,nl)
  bdirsz = fltarr(nz,nf,nl)

  close,1
  openr,1,dir+'bdirsxu.dat',/f77_unformatted
  readu,1,bdirsx
  close,1
  openr,1,dir+'bdirsyu.dat',/f77_unformatted
  readu,1,bdirsy
  close,1
  openr,1,dir+'bdirszu.dat',/f77_unformatted
  readu,1,bdirsz

; gsr

  gsrx = fltarr(nz,nf,nl)
  gsry = fltarr(nz,nf,nl)
  gsrz = fltarr(nz,nf,nl)

  close,1
  openr,1,dir+'gsrxu.dat',/f77_unformatted
  readu,1,gsrx
  close,1
  openr,1,dir+'gsryu.dat',/f77_unformatted
  readu,1,gsry
  close,1
  openr,1,dir+'gsrzu.dat',/f77_unformatted
  readu,1,gsrz

; gstheta

  gsthetax = fltarr(nz,nf,nl)
  gsthetay = fltarr(nz,nf,nl)
  gsthetaz = fltarr(nz,nf,nl)

  close,1
  openr,1,dir+'gsthetaxu.dat',/f77_unformatted
  readu,1,gsthetax
  close,1
  openr,1,dir+'gsthetayu.dat',/f77_unformatted
  readu,1,gsthetay
  close,1
  openr,1,dir+'gsthetazu.dat',/f77_unformatted
  readu,1,gsthetaz

; gsphi

  gsphix = fltarr(nz,nf,nl)
  gsphiy = fltarr(nz,nf,nl)
  gsphiz = fltarr(nz,nf,nl)

  close,1
  openr,1,dir+'gsphixu.dat',/f77_unformatted
  readu,1,gsphix
  close,1
  openr,1,dir+'gsphiyu.dat',/f77_unformatted
  readu,1,gsphiy
  close,1
  openr,1,dir+'gsphizu.dat',/f77_unformatted
  readu,1,gsphiz

; vector

  xrg   = fltarr(nz,nf,nl)
  xthg  = fltarr(nz,nf,nl)
  xphig = fltarr(nz,nf,nl)

  close,1
  openr,1,dir+'xrgu.dat',/f77_unformatted
  readu,1,xrg
  close,1
  openr,1,dir+'xthgu.dat',/f77_unformatted
  readu,1,xthg
  close,1
  openr,1,dir+'xphigu.dat',/f77_unformatted
  readu,1,xphig


; vps

  vpsnx = fltarr(nz,nf,nl)
  vpsny = fltarr(nz,nf,nl)
  vpsnz = fltarr(nz,nf,nl)

  close,1
  openr,1,dir+'vpsnxu.dat',/f77_unformatted
  readu,1,vpsnx
  close,1
  openr,1,dir+'vpsnyu.dat',/f77_unformatted
  readu,1,vpsny
  close,1
  openr,1,dir+'vpsnzu.dat',/f77_unformatted
  readu,1,vpsnz

; vhs

  vhsnx = fltarr(nz,nf,nl)
  vhsny = fltarr(nz,nf,nl)
  vhsnz = fltarr(nz,nf,nl)

  close,1
  openr,1,dir+'vhsnxu.dat',/f77_unformatted
  readu,1,vhsnx
  close,1
  openr,1,dir+'vhsnyu.dat',/f77_unformatted
  readu,1,vhsny
  close,1
  openr,1,dir+'vhsnzu.dat',/f77_unformatted
  readu,1,vhsnz

; parallel to B 

  gparb = gsrx * bdirsx + gsry * bdirsy + gsrz * bdirsz

; perpendicular to B 

  gperb = gsrx * vpsnx + gsry * vpsny + gsrz * vpsnz

; perp to B 

  gperbh = gsrx * vhsnx + gsry * vhsny + gsrz * vhsnz

  end

