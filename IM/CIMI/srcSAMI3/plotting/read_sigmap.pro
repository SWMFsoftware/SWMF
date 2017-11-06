

; pedersen conductivity data

  close,1

  sigmap   = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'sigmapu.dat',/f77_unformatted

  print,' reading sigmap '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  sigmap(*,*,0:nl-1,*) = tmpt
  sigmap(*,*,nl,*)     = sigmap(*,*,0,*)

  sigmappp               = fltarr(nz,nf+1,nl+1,nt)
  sigmappp(*,0:nf-2,*,*) = sigmap(*,0:nf-2,*,*)

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      sigmaptot = 0.
      for k = 0,nl do begin
        sigmaptot  = sigmap(i,nf-2,k,n) + sigmaptot
     endfor
     sigmappp(i,nf-1,0,n) = sigmaptot / float(nl)
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      for k = 1,nl do begin
        sigmappp(i,nf-1,k,n) = sigmappp(i,nf-1,0,n)
      endfor
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      sigmaptot = 0.
      for k = 0,nl do begin
        sigmaptot  = sigmap(i,nf-1,k,n) + sigmaptot
     endfor
     sigmappp(i,nf,0,n) = sigmaptot / float(nl)
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      for k = 1,nl do begin
        sigmappp(i,nf,k,n) = sigmappp(i,nf,0,n)
      endfor
    endfor  
  endfor 

  close,1

  delvar,tmpr,etmp

  end
