

; hall conductivity data

  close,1

  sigmah   = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'sigmahu.dat',/f77_unformatted

  print,' reading sigmah '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  sigmah(*,*,0:nl-1,*) = tmpt
  sigmah(*,*,nl,*)     = sigmah(*,*,0,*)

  sigmahpp               = fltarr(nz,nf+1,nl+1,nt)
  sigmahpp(*,0:nf-2,*,*) = sigmah(*,0:nf-2,*,*)

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      sigmahtot = 0.
      for k = 0,nl do begin
        sigmahtot  = sigmah(i,nf-2,k,n) + sigmahtot
     endfor
     sigmahpp(i,nf-1,0,n) = sigmahtot / float(nl)
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      for k = 1,nl do begin
        sigmahpp(i,nf-1,k,n) = sigmahpp(i,nf-1,0,n)
      endfor
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      sigmahtot = 0.
      for k = 0,nl do begin
        sigmahtot  = sigmah(i,nf-1,k,n) + sigmahtot
     endfor
     sigmahpp(i,nf,0,n) = sigmahtot / float(nl)
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      for k = 1,nl do begin
        sigmahpp(i,nf,k,n) = sigmahpp(i,nf,0,n)
      endfor
    endfor  
  endfor 

  close,1

  delvar,tmpr,etmp

  end
