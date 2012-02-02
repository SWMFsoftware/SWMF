
pro gitm_read_bin, file, data, time, nVars, Vars, version

  filelist = findfile(file)

  nFiles = n_elements(filelist)

  if (nFiles gt 1) then Time = dblarr(nFiles) else Time = 0.0D

  for iFile = 0, nFiles-1 do begin

      filein = filelist(iFile)

      close, 1
      openr, 1, filein, /f77

      version = 0.0D

      nLons = 0L
      nLats = 0L
      nAlts = 0L
      nVars = 0L

      readu, 1, version
      readu, 1, nLons, nLats, nAlts
      readu, 1, nVars

      Vars = strarr(nVars)
      line = bytarr(40)
      for iVars = 0, nVars-1 do begin
          readu, 1, line
          Vars(iVars) = strcompress(string(line),/remove)
      endfor

      lTime = lonarr(7)
      readu, 1, lTime

      iTime = fix(lTime(0:5))
      c_a_to_r, itime, rtime
      Time(iFile) = rTime + lTime(6)/1000.0

      if (nFiles eq 1) then begin
          Data = dblarr(nVars, nLons, nLats, nAlts)
      endif else begin
          if (iFile eq 0) then $
            Data = dblarr(nFiles, nLons, nLats, nAlts, nVars)
      endelse

      tmp = dblarr(nLons, nLats, nAlts)
      for i=0,nVars-1 do begin
          readu,1,tmp
          data(i,*,*,*) = tmp
      endfor
          
      close, 1

  endfor

end

;;filelist = findfile('3DALL*.bin')
;;file = filelist(1)
;;
;;gitm_read_bin, file, data, time, nVars, Vars, version
;;
;;contour, alog10(data(*,*,10,26)+1.0e-5), /follow, nlevels = 18, /fill
;;;contour, reform(data(*,0,0,*,18)), /follow, nlevels = 18, /fill
;;
;;end
