; READ_TEC.PRO

  close,1
  blat0r = fltarr(nx,ny,nl)
  openr,1,dir+'blat0.dat',/f77_unformatted
  print,' reading blat0'
  readu,1,blat0r

  blat0 = fltarr(nx,ny,nl+1)
  blat0(*,*,0:nl-1 ) = blat0r(*,*,*)
  blat0(*,*,nl)      = blat0r(*,*,0)

  close,1
  blon0r = fltarr(nx,ny,nl)
  openr,1,dir+'blon0.dat',/f77_unformatted
  print,' reading blon0'
  readu,1,blon0r

  blon0 = fltarr(nx,ny,nl+1)
  blon0(*,*,0:nl-1 ) = blon0r(*,*,*)
  blon0(*,*,nl)      = blon0r(*,*,0)

  close,1
  balt0r = fltarr(nx,ny,nl)
  openr,1,dir+'balt0.dat',/f77_unformatted
  print,' reading balt0'
  readu,1,balt0r

  balt0 = fltarr(nx,ny,nl+1)
  balt0(*,*,0:nl-1 ) = balt0r(*,*,*)
  balt0(*,*,nl)      = balt0r(*,*,0)

  close,1
  tecr = fltarr(nx,nl,nt)
  openr,1,dir+'tecub.dat',/f77_unformatted
  print,' reading tecb'
  readu,1,tecr

  tecB = fltarr(nx,nl+1,nt)
  tecB(*,0:nl-1,*) = tecr(*,*,*)
  tecB(*,nl,*)     = tecr(*,0,*)

  close,1
  nmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'nmf2ub.dat',/f77_unformatted
  print,' reading nmf2b'
  readu,1,nmf2r

  nmf2B = fltarr(nx,nl+1,nt)
  nmf2B(*,0:nl-1,*) = nmf2r(*,*,*)
  nmf2B(*,nl,*)     = nmf2r(*,0,*)

  close,1
  hmf2r = fltarr(nx,nl,nt)
  openr,1,dir+'hmf2ub.dat',/f77_unformatted
  print,' reading hmf2b'
  readu,1,hmf2r

  hmf2B = fltarr(nx,nl+1,nt)
  hmf2B(*,0:nl-1,*) = hmf2r(*,*,*)
  hmf2B(*,nl,*)     = hmf2r(*,0,*)

  close,1

  end
