
pro c_ymd_to_a, itime, ymd

  itime = intarr(6)

  itime(0) = fix(strmid(ymd,0,4))
  itime(1) = fix(strmid(ymd,4,2))
  itime(2) = fix(strmid(ymd,6,2))

end
