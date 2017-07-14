
pro determine_min_max, data, mini, maxi, percent = percent

  if (n_elements(percent) eq 0) then percent = 0.995 $
  else percent = float(percent)/100

  mini = min(data)
  maxi = max(data)
  r = maxi-mini
  dr = r/1000.0

  ; Determine maximum:
  n = float(n_elements(data))
  l = where(data lt maxi,c)
  p = float(c)/n
  while p gt percent do begin
     maxi = maxi - dr
     l = where(data lt maxi,c)
     p = float(c)/n
  endwhile

  ; Determine minimum:
  n = float(n_elements(data))
  l = where(data gt mini,c)
  p = float(c)/n
  while p gt percent do begin
     mini = mini + dr
     l = where(data gt mini,c)
     p = float(c)/n
  endwhile


end

