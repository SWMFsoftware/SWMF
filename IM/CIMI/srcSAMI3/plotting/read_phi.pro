
; READ_PHI.PRO

; phi

  close,1
  phi    = fltarr(nnx,nny,nt)
  phitmp = fltarr(nnx,nny)
  openr,1,dir+'phiu.dat',/f77_unformatted
  print,' reading phi'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,phitmp
    phi(*,*,i) = phitmp
  endfor

  phip    = fltarr(nf,nl+1,nt)
  phip(0:nf-2,*,*) = transpose(phi(*,0:nf-2,*),[1,0,2])
  dbang = (blat(nz-1,nf-1,0)-90.) / (blat(nz-1,nf-2,0)-90.)
  print,'dbang',dbang
  phip(nf-1,*,*)   = phip(nf-2,*,*) * dbang

  phipp              = fltarr(nf+1,nl+1,nt)
  phipp(0:nf-1,*,*)  = phip(0:nf-1,*,*)
  phipp(nf,*,*)      = 0.

  blatpp             = fltarr(nz,nf+1,nl+1)
  blatpp(*,0:nf-1,*) = blat(*,0:nf-1,*)
  blatpp(*,nf,*)     = 90.

  blonpp             = fltarr(nz,nf+1,nl+1)
  blonpp(*,0:nf-1,*) = blon(*,0:nf-1,*)
  blonpp(*,nf,*)     = blon(*,nf-1,*)


  close,1

  end







