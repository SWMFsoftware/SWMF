
; READ_MAG_COORD.PRO

; magnetic longitude

  close,1
  blonr = fltarr(nz,nf,nl)
  openr,1,dir+'blonu.dat',/f77_unformatted
  print,' reading blon'
  readu,1,blonr

  blon = fltarr(nz,nf,nl+1)
  blon(*,*,0:nl-1) = blonr
  blon(*,*,nl)   = blon(*,*,0)

; magnetic longitude

  close,1
  blonpr = fltarr(nz+1,nf+1,nl)
  openr,1,dir+'blonpu.dat',/f77_unformatted
  print,' reading blonp'
  readu,1,blonpr

  blonp = fltarr(nz+1,nf+1,nl+1)
  blonp(*,*,0:nl-1) = blonpr
  blonp(*,*,nl)   = blonp(*,*,0)

; magnetic latitude

  close,1
  blatr = fltarr(nz,nf,nl)
  openr,1,dir+'blatu.dat',/f77_unformatted
  print,' reading blat'
  readu,1,blatr

  blat = fltarr(nz,nf,nl+1)
  blat(*,*,0:nl-1) = blatr
  blat(*,*,nl)   = blat(*,*,0)

; magnetic latitude

  close,1
  blatpr = fltarr(nz+1,nf+1,nl)
  openr,1,dir+'blatpu.dat',/f77_unformatted
  print,' reading blatp'
  readu,1,blatpr

  blatp = fltarr(nz+1,nf+1,nl+1)
  blatp(*,*,0:nl-1) = blatpr
  blatp(*,*,nl)   = blatp(*,*,0)

; altitude

  close,1
  baltr = fltarr(nz,nf,nl)
  openr,1,dir+'baltu.dat',/f77_unformatted
  print,' reading balt'
  readu,1,baltr

  balt = fltarr(nz,nf,nl+1)
  balt(*,*,0:nl-1) = baltr - 6370.
  balt(*,*,nl)   = balt(*,*,0)

; altitude

  close,1
  baltpr = fltarr(nz+1,nf+1,nl)
  openr,1,dir+'baltpu.dat',/f77_unformatted
  print,' reading baltp'
  readu,1,baltpr

  baltp = fltarr(nz+1,nf+1,nl+1)
  baltp(*,*,0:nl-1) = baltpr - 6370.
  baltp(*,*,nl)   = baltp(*,*,0)

  end
