
filelist = findfile('log*.dat')

nFiles = n_elements(filelist)

nLines = 0L

for i=0,nFiles-1 do begin

    spawn, 'wc '+filelist(i), nL
    nLines = nLines+long(nL)-2

endfor

line = ''

iLine = 0L

for i=0,nFiles-1 do begin

    openr, 1, filelist(i)

    while (strpos(line,"START") lt 0) do readf,1,line
    readf,1,line

    if (i eq 0) then begin

        l = strsplit(line)
        nVars = n_elements(l)
        l = [l,strlen(line)]
        Vars = strarr(nVars)
      
        for iVar = 0, nVars-1 do Vars(iVar) = strmid(line,l(iVar),l(iVar+1)-l(iVar))

        data = fltarr(nVars,nLines)
        tmp = fltarr(nVars)

    endif

    while not eof(1) do begin

        readf,1,tmp
        data(*,iLine) = tmp
        iLine = iLine + 1L

    endwhile

    close,1

endfor

setdevice, 'log.ps','p',5

nLines = iLine
time = dblarr(nLines)

for i=0L,nLines-1 do begin

    itime = fix(data(1:6,i))
    c_a_to_r, itime, rtime
    time(i) = rtime

endfor

stime = time(0)
c_r_to_a, itime, stime
c_a_to_s, itime, strtime

t = time-time(0)
cTime = 'Simulation Time'

if (max(t) gt 7200.0) then begin
    t = t/3600.0
    cTime = cTime+' (Hours)'
    nDays = max(t)/(24.0)
endif else begin
    cTime = cTime+' (Seconds)'
    nDays = 0
endelse

plot, t, data(10,*), linestyle = 2, ytitle = 'Temperature (K)', $
  xtitle = cTime, xstyle = 1, pos = [0.1,0.3,0.9,0.7]
oplot, t, data(9,*), linestyle = 2
oplot, t, data(11,*), thick = 3

xyouts, 0.9, 0.705, /norm, 'Start Time : '+strtime, alignment = 1

if (nDays gt 1) then begin
   for i=1, nDays do begin
      oplot, [i,i]*24.0, [0, max(t)*10.0], linestyle = 1
   endfor
endif

closedevice

end


