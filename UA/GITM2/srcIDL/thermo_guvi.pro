filelist=findfile('*sav')
nfiles=n_elements(filelist)
print,'nfiles=',nfiles

restore,filelist(0)

outfile=''
outfile='GUVI_'+strmid(saved_data[0].year,2,2)+ $
         string(saved_data[0].month,format='(I2)')+ $
         string(saved_data[0].day,format='(I2)')+'.dat'
print,outfile

norbits=n_elements(saved_data)

openw,1,'temp.dat', /append

for iorbit=norbits-1,0,-1 do begin
print,'iorbit=',iorbit
ntotal=n_elements(saved_data[iorbit].UT)
done=0
istart=0
iend=1

timearray=intarr(3,3468)
timearray(0,*)=floor(saved_data[iorbit].UT(*))
timearray(1,*)=floor((saved_data[iorbit].UT(*)-timearray(0,*))*60.)
timearray(2,*)=((saved_data[iorbit].UT(*)-timearray(0,*))*60.-timearray(1,*))*60.


while not done do begin
fmt='(7I5,4f9.2)'
if saved_data[iorbit].UT(iend) lt saved_data[iorbit].UT(iend-1) then begin
icenter = (istart+iend)/2
if saved_data[iorbit].ut(icenter) gt 0 then $
printf,1,fix(saved_data[iorbit].year),fix(saved_data[iorbit].month),fix(saved_data[iorbit].day),$
         timearray(0:2,icenter),0,saved_data[iorbit].lon(icenter),saved_data[iorbit].lat(icenter),$
         0.0,saved_data[iorbit].on2(icenter),format=fmt

istart=iend+1
iend=istart
endif
iend=iend+1
if iend eq ntotal then done=1
if saved_data[iorbit].ut(iend) eq 0 then done=1
endwhile
endfor
close,1

data1=intarr(7,3468)
data2=fltarr(4,3468)
tmp1=intarr(7)
tmp2=fltarr(4)
i=0

openr,2,'temp.dat'
while not eof(2) do begin
readf,2,tmp1,tmp2,format=fmt
data1(*,i)=tmp1(*)
data2(*,i)=tmp2(*)
i=i+1
endwhile
close,2
print,'number of line',i

openw,3,outfile
printf,3,'#start'
for j=i-1,0,-1 do begin
printf,3, data1(*,j),data2(*,j),format=fmt
endfor
close,3

data_2D=reform(data2,4,4,867)

contour,data_2D(3,*,*),data_2D(0,*,*),data_2D(1,*,*),/fill,levels=findgen(21)/20

end


;openr,1,'ON2_2003_302m.sav'
;
;ON2_GRID=fltarr(204,102)
;TMP1 = fltarr(204)
;TMP2 = {ON2:fltarr(3468),LAT:fltarr(3468),LON:fltarr(3468),UT:fltarr(3468),$
;       YEAR:'2003',DOY:'302',ORBIT:'10227',points:2089,month:10,day:29,D_lat:1.76471,D_lon:1.76471,sza:fltarr(3468)}
;saved_data=replicate(tmp2,15)
;
;for i=0,101 do begin
;print,i
;readf,1,tmp1
;ON2_GRID(*,i)=tmp1
;endfor
;
; i = 0
;;readf,1,tmp1
; while not EOF(1)do begin
;print,i
;readf,1,tmp2
;saved_data(i)=tmp2
;i=i+1
;ENDWHILE
;print,'read file',i
;CLOSE,1         
;



