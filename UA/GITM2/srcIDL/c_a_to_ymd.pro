
pro c_a_to_ymd, iTime, ymd, twoyear = twoyear

  iTimeIn = iTime

  c_a_to_r, iTimeIn, rTime
  c_r_to_a, iTimeIn, rTime

  if (n_elements(twoyear) gt 0) then chop = 2 else chop = 4
  sYear  = chopr(tostr(iTimeIn(0)),chop)
  sMonth = chopr('00'+tostr(iTimeIn(1)),2)
  sDay   = chopr('00'+tostr(iTimeIn(2)),2)

  ymd    = sYear+sMonth+sDay

end

