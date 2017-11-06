

; electron density data

  close,1

  dene   = fltarr(nz,nf,nl,nt)
  tmp0   = fltarr(nz,nf,nl)

  openr,1,dir+'deneu.dat',/f77_unformatted

  print,' reading dene '

  for i = 0,nt-1 do begin
    print,i 
    readu,1,tmp0
    dene(*,*,*,i) = tmp0
  endfor

;  dene(*,*,0:nl-1,*) = tmpt
;  dene(*,*,nl,*)     = dene(*,*,0,*)

;  denepp               = fltarr(nz,nf+1,nl+1,nt)
;  denepp(*,0:nf-2,*,*) = dene(*,0:nf-2,*,*)

;  for n = 0,nt-1 do begin
;    for i = 0,nz-1 do begin
;      denetot = 0.
;      for k = 0,nl do begin
;        denetot  = dene(i,nf-2,k,n) + denetot
;     endfor
;     denepp(i,nf-1,0,n) = denetot / float(nl)
;    endfor  
;  endfor 

;  for n = 0,nt-1 do begin
;    for i = 0,nz-1 do begin
;      for k = 1,nl do begin
;        denepp(i,nf-1,k,n) = denepp(i,nf-1,0,n)
;      endfor
;    endfor  
;  endfor 

;  for n = 0,nt-1 do begin
;    for i = 0,nz-1 do begin
;      denetot = 0.
;      for k = 0,nl do begin
;        denetot  = dene(i,nf-1,k,n) + denetot
;     endfor
;     denepp(i,nf,0,n) = denetot / float(nl)
;    endfor  
;  endfor 

;  for n = 0,nt-1 do begin
;    for i = 0,nz-1 do begin
;      for k = 1,nl do begin
;        denepp(i,nf,k,n) = denepp(i,nf,0,n)
;      endfor
;    endfor  
;  endfor 

  close,1

  delvar,tmpr,etmp

  end
