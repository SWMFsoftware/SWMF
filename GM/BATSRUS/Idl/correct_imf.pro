pro correct_imf,wIn,xIn,inputfile,outputfile

; wIn contains the upstream data with 15 columns:
;    yr mo dy hr mn sc ms bx by bz ux uy uz rho T
; xIn should be set to the X position of the satellite 
     relative to the inflow boundary in units of km !
; inputfile is the name of the original IMF file 
; outputfile is the name of the corrected IMF file

w = wIn
x = xIn

; Check arrays
siz = size(w)
if siz(0) ne 2 then begin
   print,'ERROR: w array should be 2D! Size(w)=',siz
   retall
endif
nPoint = siz(1)
nVar   = siz(2)

if nVar ne 15 then begin
   print,'ERROR: nVar (2nd dimension of w) is not 15! nVar=',nVar
   retall
endif

if n_elements(x) ne nPoint then begin
    print,'ERROR: w and x arrays have different sizes:',$
      nPoint,n_elements(x)
    retall
endif

; Calculate the time in seconds measured from the beginning of the month
iyr=0 & imo=1 & idy = 2 & ihr = 3 & imn = 4 & isc = 5 & ims = 6

Time = $
  w(*,idy)*86400. + $
  w(*,ihr)* 3600. + $
  w(*,imn)*   60. + $
  w(*,isc)        + $
  w(*,ims)*    0.001

; Calculate the time delay
iux = 10
TimeDelay = -x/w(*,iux)
NewTime   = Time+TimeDelay


; Smooth out rarefaction waves
ibx = 7
for iVar=ibx,nVar-1 do begin
    i0   = 0
    Var0 = w(0,iVar)
    for i1=1,nPoint-1 do begin
        if w(i1,iVar) ne Var0 then begin
                                ; For rarefaction waves the negative
                                ; velocity is increasing
                                ; Maybe a more sophisticated shock
                                ; recognition is needde

            if w(i1,iux) gt w(i1-1,iux) then begin
                for i=i0+1,i1-1 do begin
                    w(i,iVar)= ((Time(i)-Time(i0))*w(i1,iVar) + $
                                (Time(i1)-Time(i))*w(i0,iVar)   $
                               )/(Time(i1)-Time(i0))
                endfor
            endif
            i0   = i1
            Var0 = w(i1,iVar)
        endif
    endfor
endfor

; Take shock waves into account
for i=1,nPoint-1 do begin
                                ; Set time=-1 for previous 
                                ; entries with longer arrival time    
    index = where(NewTime(0:i-1) gt NewTime(i))
    if index(0) ge 0 then NewTime(index) = -1
endfor

; Print out corrected IMF file
close,1
openw,1,outputfile
printf,1,'Corrected IMF based on ',inputfile
printf,1,'yr mo dy hr mn sc ms bx by bz ux uy uz rho T'
for i=0,nPoint-1 do begin
    ; Skip entries marked with negative times
    if(NewTime(i) ge 0)then begin
        ; Calculate integer times
        t = NewTime(i)
        day  = fix(t/86400.) & t = t - day *86400.
        hour = fix(t/3600.)  & t = t - hour*3600.
        min  = fix(t/60.)    & t = t - min*60.
        sec  = fix(t)        & t = t - sec
        msc  = fix(1000.*t)
        printf,1,w(i,iyr),w(i,imo),day,hour,min,sec,msc,w(i,ibx:nVar-1),$
          format='(i5,5i3,i4,7f11.2,f13.2)'
    endif
endfor
close,1

end
