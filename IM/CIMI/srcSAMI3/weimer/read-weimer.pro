; READ.PRO; set the color table
  
  device,true=24
  device,de=0
  device,retain=2
  loadct,41

  nz = 204
;  nf = 124
  nf =  120
  nl = 96
  nlat = nf+1
  nlon = nl+1

  nt   = 600
; if reading only one day, nskip should be equal to nweimer.rst 
; at the end of the previous day.
  nskip = 0

  nzp1 = nz+1
  nfp1 = nf+1
  nfm1 = nf-1
  nlp1 = nl+1

;  blonp = fltarr(nzp1,nfp1,nlp1)
;  blatp = fltarr(nzp1,nfp1,nlp1)

;  close,1,2
;  openr,1,dir+'blonpu.dat',/f77_unformatted
;  openr,2,dir+'blatpu.dat',/f77_unformatted

;  readu,1,blonp
;  readu,2,blatp
;  close,1,2

  timew = fltarr(5,nt)

  epot = fltarr(nlat,nlon,nt)
  tmp  = dblarr(nlat,nlon)
  hrut = fltarr(nt)
  fac  = fltarr(nlat,nlon,nt)

;  openr,1,dir+'fac_weimer.inp',/f77_unformatted
;  for i = 0,nt-1 do begin
;    readu,1,x
;    print,'i,x',i,x
;    readu,1,tmp

;    fac(*,*,i) = tmp
;  endfor

  dir='./'

  close,1
  x=1.1d
  openr,1,dir+'phi_weimer.inp',/f77_unformatted

  for i = 0,nskip-1 do begin
    readu,1,x
    readu,1,tmp
  endfor
  for i = 0,nt-1 do begin
    readu,1,x
    hrut(i)     = float(x)
    readu,1,tmp
    epot(*,*,i) = float(tmp)
  endfor
  close,1

; get day and year from weimer_input.inp
  close,1
  openr,1,dir+'weimer_input.inp' 
  a = ' '
  readf,1,a
  readf,1,a
  readf,1,a
  readf,1,a
  n1=1
  nyear=0
  nday=0
  readf,1,n1,nyear,nday

  print,'n1,nyear,nday',n1,nyear,nday

;  construct timew, corresponds to usual time file
  for k=0,nt-1 do begin
    timew(0,k) = k+1
    timew(1,k) = fix(hrut(k))
    timew(2,k) = fix(60.*(hrut(k) - fix(hrut(k))))
    timew(3,k) = 3600.0*(hrut(k) - timew(1,k) - timew(2,k)/60.0)
    timew(4,k) = hrut(k) 
    timew(1,k) = timew(1,k) MOD 24
  endfor

;  for i=0,6 do begin
;    for k=0,nt-1 do begin
;      if (timew(4,k) lt timew(4,k-1)) timew(4,k) = timew(4,k) + 24.0
;    endfor
;  endfor    

; phi from SAMI3 has dimensions        nlp1,nfm1,nt
; blon,blat from sami3 have dimensions nz  ,nf  ,nlp1
; map-phi.pro uses phi(*,0:nf-2,ntm) transposed
;                  blon(nz-1,*,0:nf-2)
;                  blat(nz-1,*,0:nf-2)
;???? does blonp(nz-1,0:nf-2,0) line up with blon(nz-1,0:nf-2,0)?

;  epot = fltarr(nfp1,nlp1,nt)
  phiw = fltarr(nlp1,nfm1,nt)
  for k = 0,nt-1 do begin
    for j = 0,nf-2 do begin
      for l = 0,nl do begin
        phiw(l,j,k) = epot(j,l,k)
      endfor
    endfor
  endfor

  end
