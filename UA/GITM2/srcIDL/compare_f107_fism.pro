
close,2
openw,2,'waves.dat'

f107file = 'srcData/f107.txt'
f107_read,f107file, f107time, f107data

dt = 81.0*86400.0
smooth_1d, f107time, f107data, dt, f107aData

filelist = findfile('srcData/FISM/fismflux_daily_200*.dat')
nFiles = n_elements(filelist)
printf,2, nFiles
for i=0,nFiles-1 do begin
   fismfile = filelist(i)
   printf,2,fismfile
   fism_read, fismfile, fismtime0, fismdata0
   if (i eq 0) then begin
      fismtime = fismtime0
      fismdata = fismdata0
   endif else begin
      fismtime = [fismtime,fismtime0]
      fismdata = [fismdata,fismdata0]      
   endelse
endfor

mintime = min(fismtime)
maxtime = max(fismtime)

l = where(f107time ge mintime and f107time le maxtime)

f107 = f107data(l)
f107a = f107adata(l)

nWaves = n_elements(fismdata(0,*))

printf,2,nWaves
for iWave = 0,nWaves-1 do begin

   fism = reform(fismdata(*,iWave))

   plot, f107, fism, ystyle = 1, xstyle = 1, psym = 4

   fit_fism_f107, fism, f107, f107a, slope, slopea, sloped, intercept, p, pa

   printf,2, iWave, slope, slopea, sloped, intercept, p, pa
   
endfor

close,2

;plot, f107data(l)-f107adata(l), fismdata(*,0), ystyle = 1, xstyle = 1, psym = 4


end
