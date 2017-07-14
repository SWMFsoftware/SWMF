pro time_axis, stime, etime, s_time_range, e_time_range, 	$
	xtickname, xtitle, xtickvalue, xminor, xtickn

  times = intarr(1,2,6)

  if n_elements(stime) eq 6 then begin
    times(0,0,*) = stime
    times(0,1,*) = etime
  endif else begin
    c_r_to_a, itime, stime
    times(0,0,*) = itime
    c_r_to_a, itime, etime
    times(0,1,*) = itime
  endelse

  if times(0,0,0) gt 1900 then times(0,0,0) = times(0,0,0) - 1900
  if times(0,0,0) ge 100 then times(0,0,0) = times(0,0,0) - 100
  if times(0,1,0) gt 1900 then times(0,1,0) = times(0,1,0) - 1900
  if times(0,1,0) ge 100 then times(0,1,0) = times(0,1,0) - 100

  placement = intarr(1,1)*0
  curvar = 0

  c_a_to_r, times(0,0,*), bt

  basetime = dblarr(1)
  basetime(0) = bt

  moncheck = "JanFebMarAprMayJunJulAugSepOctNovDec"

  compute_axis, times, placement, s_time_range, e_time_range,	$
	curvar, xtickname, xtitle, basetime, moncheck,		$
	xtickn, xtickvalue, xminor

  return

end