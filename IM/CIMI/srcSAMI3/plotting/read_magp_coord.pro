
; READ_MAGP_COORD.PRO

; magnetic longitude

  close,1
  blonp = fltarr(nz+1,nf+1,nl)
  openr,1,dir+'blonpu.dat',/f77_unformatted
  print,' reading blonp'
  readu,1,blonp

; magnetic latitude

  close,1
  blatp = fltarr(nz+1,nf+1,nl)
  openr,1,dir+'blatpu.dat',/f77_unformatted
  print,' reading blatp'
  readu,1,blatp

; altitude

  close,1
  baltp = fltarr(nz+1,nf+1,nl)
  openr,1,dir+'baltpu.dat',/f77_unformatted
  print,' reading baltp'
  readu,1,baltp

  end
