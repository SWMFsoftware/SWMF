function jday, year, mon, day
  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
  if year mod 4 eq 0 then dayofmon(1) = dayofmon(1) + 1
  julianday = 0
  for i=0,mon-2 do julianday = julianday + dayofmon(i)
  julianday = julianday + day

  return, julianday
end
