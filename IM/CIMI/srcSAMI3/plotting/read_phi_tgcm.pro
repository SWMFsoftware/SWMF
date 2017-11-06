
; read_phi_tgcm

  thesize  = size(phi_tgcm)
  jlon     = thesize(1)

  for i = 0,jlon-1 do begin
    if lons(i) lt 0  then lons(i) = lons(i) + 360.
  endfor

  end
