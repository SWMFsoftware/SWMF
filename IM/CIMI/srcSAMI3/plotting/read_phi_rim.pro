; READ_PHI.PRO

; set the color table                                                                     

  device,true=24
  device,de=0
  device,retain=2
  loadct,41

  ilon = 80
  jlat = 48

  nnx  = ilon + 1
  nny  = jlat - 1

  ylonp       = fltarr(nnx)
  ylatp       = fltarr(nny)
  phi         = fltarr(nnx,nny)

  close,1
  openr,1,dir+'lonlatphiu.dat',/f77_unformatted
  readu,1,ylonp,ylatp,phi
  close,1

; fold grid/phi to do both NS hemispheres

  ylatpp      = fltarr(2*nny)
  phip        = fltarr(nnx,2*nny)

   rim1_tmp    = fltarr(nnx,nny)
   rim1p_sami  = fltarr(nnx,2*nny)
   rim1_tmp(0:nnx-2,0:nny-1) = rim1_sami(0:nnx-2,0:nny-1)
   rim1_tmp(nnx-1,*)         = rim1_tmp(0,*)

   rim2_tmp    = fltarr(nnx,nny)
   rim2p_sami  = fltarr(nnx,2*nny)
   rim2_tmp(0:nnx-2,0:nny-1) = rim2_sami(0:nnx-2,0:nny-1)
   rim2_tmp(nnx-1,*)         = rim2_tmp(0,*)

   zigm11_tmp    = fltarr(nnx,nny)
   zigm11p_sami  = fltarr(nnx,2*nny)
   zigm11_tmp(0:nnx-2,0:nny-1) = zigm11_sami(0:nnx-2,0:nny-1)
   zigm11_tmp(nnx-1,*)         = zigm11_tmp(0,*)

   zigm22_tmp    = fltarr(nnx,nny)
   zigm22p_sami  = fltarr(nnx,2*nny)
   zigm22_tmp(0:nnx-2,0:nny-1) = zigm22_sami(0:nnx-2,0:nny-1)
   zigm22_tmp(nnx-1,*)         = zigm22_tmp(0,*)

   zigm2_tmp    = fltarr(nnx,nny)
   zigm2p_sami  = fltarr(nnx,2*nny)
   zigm2_tmp(0:nnx-2,0:nny-1) = zigm2_sami(0:nnx-2,0:nny-1)
   zigm2_tmp(nnx-1,*)         = zigm2_tmp(0,*)


  for j = nny,2*nny-1 do begin
    jj        = j - nny 
    ylatpp(j) = ylatp(jj)
  endfor

  for j = 0,nny-1 do begin
    jj        = nny-1-j
    ylatpp(j) = -ylatp(jj)
  endfor

  for j = 0,2*nny-1 do print,j,ylatpp(j)

  for j = nny,2*nny-1 do begin
    for i = 0,nnx-1 do begin
      jj        = j - nny 
      phip(i,j) = phi(i,jj)
      rim1p_sami(i,j) = rim1_tmp(i,jj)
      rim2p_sami(i,j) = rim2_tmp(i,jj)
      zigm11p_sami(i,j) = zigm11_tmp(i,jj)
      zigm22p_sami(i,j) = zigm22_tmp(i,jj)
      zigm2p_sami(i,j)  = zigm2_tmp(i,jj)
    endfor
  endfor

  for j = 0,nny-1 do begin
    for i = 0,nnx-1 do begin
      jj        = nny-1-j
      phip(i,j) = phi(i,jj)
      rim1p_sami(i,j) = rim1_tmp(i,jj)
      rim2p_sami(i,j) = rim2_tmp(i,jj)
      zigm11p_sami(i,j) = zigm11_tmp(i,jj)
      zigm22p_sami(i,j) = zigm22_tmp(i,jj)
      zigm2p_sami(i,j)  = zigm2_tmp(i,jj)
    endfor
  endfor

  end
