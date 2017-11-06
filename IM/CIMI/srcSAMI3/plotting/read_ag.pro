; READ_AG.PRO

;  nx = 200
;  ny = 200

  close,1
  glat0r = fltarr(nx,ny,nl)
  openr,1,dir+'glat0.dat',/f77_unformatted
  print,' reading glat0'
  readu,1,glat0r

  glat0 = fltarr(nx,ny,nl+1)
  glat0(*,*,0:nl-1 ) = glat0r(*,*,*)
  glat0(*,*,nl)      = glat0r(*,*,0)

  close,1
  glon0r = fltarr(nx,ny,nl)
  openr,1,dir+'glon0.dat',/f77_unformatted
  print,' reading glon0'
  readu,1,glon0r

  glon0 = fltarr(nx,ny,nl+1)
  glon0(*,*,0:nl-1 ) = glon0r(*,*,*)
  glon0(*,*,nl)      = glon0r(*,*,0)

  close,1
  zalt0r = fltarr(nx,ny,nl)
  openr,1,dir+'zalt0.dat',/f77_unformatted
  print,' reading zalt0'
  readu,1,zalt0r

  zalt0 = fltarr(nx,ny,nl+1)
  zalt0(*,*,0:nl-1 ) = zalt0r(*,*,*)
  zalt0(*,*,nl)      = zalt0r(*,*,0)

  close,1
  agr = fltarr(nx,nl,nt)
  openr,1,dir+'agu.dat',/f77_unformatted
  print,' reading ag'
  readu,1,agr

  ag = fltarr(nx,nl+1,nt)
  ag(*,0:nl-1,*) = agr(*,*,*)
  ag(*,nl,*)     = agr(*,0,*)

  close,1
  nmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'nmf2u.dat',/f77_unformatted
  print,' reading nmf2'
  readu,1,nmf2r

  nmf2 = fltarr(nx,nl+1,nt)
  nmf2(*,0:nl-1,*) = nmf2r(*,*,*)
  nmf2(*,nl,*)     = nmf2r(*,0,*)

  close,1
  hmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'hmf2u.dat',/f77_unformatted
  print,' reading hmf2'
  readu,1,hmf2r

  hmf2 = fltarr(nx,nl+1,nt)
  hmf2(*,0:nl-1,*) = hmf2r(*,*,*)
  hmf2(*,nl,*)     = hmf2r(*,0,*)

  close,1

  end
