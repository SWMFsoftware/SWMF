
; READ_GEO_COORD.PRO

; magnetic longitude

  close,1
  glonr = fltarr(nz,nf,nl)
  openr,1,dir+'glonu.dat',/f77_unformatted
  print,' reading glon'
  readu,1,glonr

  glon = fltarr(nz,nf,nl+1)
  glon(*,*,0:nl-1) = glonr
  glon(*,*,nl)   = glon(*,*,0)

; magnetic latitude

  close,1
  glatr = fltarr(nz,nf,nl)
  openr,1,dir+'glatu.dat',/f77_unformatted
  print,' reading glat'
  readu,1,glatr

  glat = fltarr(nz,nf,nl+1)
  glat(*,*,0:nl-1) = glatr
  glat(*,*,nl)   = glat(*,*,0)

; altitude

  close,1
  zaltr = fltarr(nz,nf,nl)
  openr,1,dir+'zaltu.dat',/f77_unformatted
  print,' reading zalt'
  readu,1,zaltr

  zalt = fltarr(nz,nf,nl+1)
  zalt(*,*,0:nl-1) = zaltr
  zalt(*,*,nl)   = zalt(*,*,0)

  end
