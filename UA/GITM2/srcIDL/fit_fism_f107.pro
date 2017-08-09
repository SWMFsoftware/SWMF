
pro fit_fism_f107, fism, f107, f107a, slope, slopea, sloped, intercept, p, pa
  
  p = 1.0
  pa = 1.0
  
  dx = max(f107)-min(f107)
  dy = max(fism)-min(fism)
  print, dy, dx
  slope = dy/dx
  intercept = min(fism) - slope * min(f107)

  slopea = slope
  sloped = slope/10.0

  fit = f107 * slope + f107a * slopea + (f107a-f107) * sloped + intercept

  p_save  = 1.0
  pa_save = 1.0

;  error_save = abs(mean(fit-fism)) ;
  mini = min(fism)*0.95
  error_save = sqrt(mean(((fit-fism)/(fism-mini))^2))
  print, error_save
  slope_save = slope
  slopea_save = slopea
  sloped_save = sloped
  intercept_save = intercept
  print, -1, slope, slopea, sloped, intercept, error_save

  plot, fit, fism, psym = 4, xstyle = 1, ystyle = 1
  oplot, [0,max(fism)]*2, [0,max(fism)]*2

  i = 0L
  iLast = 0L

  while (i-iLast lt 10000 and i lt 500000) do begin

     fac = 25.0
     if (i gt 1000) then fac = 100.0
     if (i gt 10000) then fac = 500.0
     if (i gt 50000) then fac = 1000.0
     
     slope      = slope_save     + randomn(s,1)*slope_save/fac
     slope = slope(0)
     slopea     = slopea_save    + randomn(s,1)*slopea_save/fac
     slopea = slopea(0)
     sloped     = sloped_save    + randomn(s,1)*sloped_save/fac
     sloped = sloped(0)
     intercept  = intercept_save + randomn(s,1)*intercept_save/fac
     intercept = intercept(0)

     p  = p_save + randomn(s,1)/100.0
     p = p(0)
     pa = pa_save + randomn(s,1)/100.0
     pa = pa(0)
     
     fit = (f107^pa) * slope + (f107a^pa) * slopea + (f107a-f107)*sloped + intercept
     ;     error = abs(mean(fit-fism))
     error = sqrt(mean(((fit-fism)/(fism-mini))^2))

;     print, i, slope, slopea, intercept, error
     
     
     if (error lt error_save) then begin
        slope_save = slope
        slopea_save = slopea
        sloped_save = sloped
        intercept_save = intercept
        p_save = p
        pa_save = pa
        error_save = error
        print, i, slope, slopea, sloped, intercept, p, pa, error
        plot, fit, fism, psym = 4, xstyle = 1, ystyle = 1
        oplot, [0,max(fism)]*2, [0,max(fism)]*2
;        wait, 0.01
        iLast = i
     endif

     i++
     
  endwhile
  
  
end
