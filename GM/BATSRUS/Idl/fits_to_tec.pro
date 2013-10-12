;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro fits_to_tec, FileIn, DataName, silent=silent

; Purpose:
;  Read fits file and write header to a file with .H extension
;  and the data into a Tecplot file.
;
; Usage:
;   fits_to_tec [,FileName] [,DataName] [,/silent]
;
; DataName is a string which contains the type of data of the image
; it is used as the data type header of the Tecplot file
; for example: DataName='Br[G]' or 'U[km/s]'
; Add /silent to suppress verbose information.

if n_elements(FileIn) eq 0 then $
  FileIn = ask('Enter input file name:','')
if n_elements(DataName) eq 0 then $
  DataName = ask('Enter data name to be written into Tecplot file:','Br[G]')

temp=strsplit(filein,'.',/extract)
FileBase  = temp(0)
if n_elements(temp) eq 1 then FileFits  = FileBase+'.fits' else FileFits=FileIn
FileHeader= FileBase+'.H'
FileTec=FileBase+'_tec.dat'

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

if not keyword_set(silent) then begin
    print,''
    print,'Writing TecPlot file ',FileTec
    print,''
endif

openw,lun,FileTec,/get_lun
printf,lun,' TITLE="',FileFits,'"'
printf,lun,'VARIABLES = "',DataName,'"'
printf,lun,'ZONE T="',filein,'", I= ',Nx,' J= ',Ny,' , K=1, F=POINT'

for i=0L,Ny-1 do for j=0L,Nx-1 do $
  printf,lun, format = '(1e14.6)',Data(j,i)

free_lun, lun

if not keyword_set(silent) then print,'Conversion done'

end
