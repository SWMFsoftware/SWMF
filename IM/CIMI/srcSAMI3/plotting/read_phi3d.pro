
; READ_PHI3D.PRO

; phi3d

  close,1
  phi3d    = fltarr(nz,nnx,nny,nt)
  phitmp = fltarr(nz,nnx,nny)
  openr,1,dir+'phi3du.dat',/f77_unformatted
  print,' reading phi3d'
  for i = 0,nt-1 do begin
    print,i 
    readu,1,phitmp
    phi3d(*,*,*,i) = phitmp
  endfor

  phip3d    = fltarr(nz,nf,nl+1,nt)
  phip3d(*,0:nf-2,*,*) = transpose(phi3d(*,*,0:nf-2,*),[0,2,1,3])
  dbang = (blat(nz-1,nf-1,0)-90.) / (blat(nz-1,nf-2,0)-90.)
  print,'dbang',dbang
  phip3d(*,nf-1,*,*)   = phip3d(*,nf-2,*,*) * dbang

;  phipp              = fltarr(nf+1,nl+1,nt)
;  phipp(0:nf-1,*,*)  = phip(0:nf-1,*,*)
;  phipp(nf,*,*)      = 0.

;  blatpp             = fltarr(nz,nf+1,nl+1)
;  blatpp(*,0:nf-1,*) = blat(*,0:nf-1,*)
;  blatpp(*,nf,*)     = 90.

;  blonpp             = fltarr(nz,nf+1,nl+1)
;  blonpp(*,0:nf-1,*) = blon(*,0:nf-1,*)
;  blonpp(*,nf,*)     = blon(*,nf-1,*)


  close,1

  end







