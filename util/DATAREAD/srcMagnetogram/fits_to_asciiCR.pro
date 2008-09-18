pro fits_to_ascii, FileIn, DataName, silent=silent,CR=CR

; Purpose:
;  Read fits file and write header to a file with .H extension
;  and the data into a ASCII file.
;
; Usage:
;   fits_to_tec [,FileName] [,DataName] [,/silent]
;
; DataName is a string which contains the type of data of the image
; it is used as the data type header of the Tecplot file
; for example: DataName='Br[G]' or 'U[km/s]'
; Add /silent to suppress verbose information.

;if n_elements(FileIn) eq 0 then begin 
;    FileIn = 'fitsfile.fits'
;endif

;nMax=180
;CR=1900

;read,nMax,prompt='enter order of harmonics (nMax:)' 
;nMax=strtrim(nMax,2)
nMax=90
;read,CR,prompt='enter Carrington Rotation number:' 
;CR=strtrim(CR,2)
;CR=2025
;
FileFits  = 'fitsfile.fits'
FileHeader='fitsfile.H'
FileDat='fitsfile.dat'
FileTec='fitsfile_tec.dat' 
DataName='Br [G]'  

Data = readfits(FileFits, ImHeader, silent=silent)

if not keyword_set(silent) then begin
    print,''
    print,'Writing header file ',FileHeader
    print,''
endif

openw,lun,FileHeader,/get_lun
printf,lun,ImHeader
free_lun, lun

; Get image dimensions
s=size(Data)
Nx=s(1)
Ny=s(2)

; Removing missing data by multiply B by sin(lat)^8
for i=0L,Ny-1 do begin
    theta=!PI*float(i)/float(Ny)
    for j=0L,Nx-1 do begin
        if(abs(Data(i*Nx+j)) ge 5000.0) then $
          Data(i*Nx+j)=Data(i*Nx+j)*sin(theta)^8
    endfor
endfor

if not keyword_set(silent) then begin
    print,''
    print,'Writing TecPlot file ',FileDat
    print,''
endif

openw,lun,FileDat,/get_lun
printf,lun,'#CR'
printf,lun,CR
printf,lun,'#nMax'
printf,lun,nMax
printf,lun,'#ARRAYSIZE'
printf,lun,strtrim(Nx,2)
printf,lun,strtrim(Ny,2)
printf,lun,'#START'

for i=0L,Ny-1 do begin
    for j=0L,Nx-1 do begin
        if(abs(Data(i*Nx+j)) gt 1900.0)then $
          Data(i*Nx+j)=abs(Data(i*Nx+j))*Data(i*Nx+j)/abs(Data(i*Nx+j)+1e-3)
        printf,lun, format = '(1e14.6)',Data(i*Nx+j)
    endfor
endfor

free_lun, lun

openw,lun,FileTec,/get_lun
printf,lun,' TITLE="',FileFits,'"'
printf,lun,'VARIABLES = "',DataName,'"'
printf,lun,'ZONE T="',FileTec,'", I= ',Nx,' J= ',Ny,' , K=1, F=POINT'

for i=0L,Ny-1 do for j=0L,Nx-1 do $
  printf,lun, format = '(1e14.6)',Data(j,i)

free_lun, lun

if not keyword_set(silent) then print,'Conversion done'

end

