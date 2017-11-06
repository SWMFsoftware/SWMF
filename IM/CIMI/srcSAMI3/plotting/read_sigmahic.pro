

; (ion) hall conductivity data 

  close,1

  sigmahic   = fltarr(nz,nf,nl+1,nt)
  tmpt   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'sigmahicu.dat',/f77_unformatted

  print,' reading sigmahic '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    tmpt(*,*,*,i) = tmp0
  endfor

  sigmahic(*,*,0:nl-1,*) = tmpt
  sigmahic(*,*,nl,*)     = sigmahic(*,*,0,*)

  sigmahicpp               = fltarr(nz,nf+1,nl+1,nt)
  sigmahicpp(*,0:nf-2,*,*) = sigmahic(*,0:nf-2,*,*)

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      sigmahictot = 0.
      for k = 0,nl do begin
        sigmahictot  = sigmahic(i,nf-2,k,n) + sigmahictot
     endfor
     sigmahicpp(i,nf-1,0,n) = sigmahictot / float(nl)
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      for k = 1,nl do begin
        sigmahicpp(i,nf-1,k,n) = sigmahicpp(i,nf-1,0,n)
      endfor
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      sigmahictot = 0.
      for k = 0,nl do begin
        sigmahictot  = sigmahic(i,nf-1,k,n) + sigmahictot
     endfor
     sigmahicpp(i,nf,0,n) = sigmahictot / float(nl)
    endfor  
  endfor 

  for n = 0,nt-1 do begin
    for i = 0,nz-1 do begin
      for k = 1,nl do begin
        sigmahicpp(i,nf,k,n) = sigmahicpp(i,nf,0,n)
      endfor
    endfor  
  endfor 

  close,1

  delvar,tmpr,etmp

  end
