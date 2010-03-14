;==========================================================================
; Written by Rona Oran , 14 Mar 2010
;
; Purpose:
;
;  1. Read 3 EIT EUV fits files (corresponding to 3 wavelengths) and merge
;  into one data set
; 
;  2. Resize image to 512 x 512 (same as synthesized images) 
;  
;  3. Write data into a .dat ASCII file, in the same
;  format that synthesized LOS EUV are outputted by BATSRUS (LC component).
;
; Usage:
;   .r euv_fits_to_ascii
;
;   you will be prompted to enter three file names, each should be a
;   fits file of an EIT image of a different wavelength. Otherwise
;   routine will be aborted with a warning.
;   
;   noresize=1
;   .r euv_fits_to_ascii
;
;   This option allows to output the combined EUV data without
;   resizing the image.  This is mostly intended for testing purposes.
;
; Output:
;
;   A file called los_euv_xxxx.dat (in TECplot format) will be
;   created ('xxxx' stands for the observation date as read from fits file)
;   
;=========================================================================
print, ''
print, 'This routine creates a single TEC(.dat) file from 3 EUV fits files.'
print, 'Use sxr_fits_to_ascii for X-Ray images.'
print, ''

; prompt user for filenames
nfile = 3
filenames=strarr(nfile)
for ifile=0,nfile-1 do begin
   print, 'Enter input filename ',ifile
   askstr,'filename ',filename,1
   filenames(ifile)=filename
   filename=''
endfor

; Assign arrays for input
npix=fltarr(ifile)
WaveLengths = strarr(nfile)
ObsDates = strarr(nfile)
ObsFullDates = strarr(nfile)

; read fits files, extract some header info 
for ifile=0,nfile-1 do begin
   Data = readfits(filenames(ifile), Header,/silent)
   print, ''
   print, 'Processing fits file ', filenames(ifile)
  
   ; Get image dimensions
   s=size(Data)
   Nx=s(1)
   Ny=s(2)
   npix(ifile) = Nx*Ny
   if Nx ne Ny then begin
      print,'WARNING!! Fits image must be square.'
      print, 'Aborting...'
      retall
   endif
   ; Get EUV wavelength
   tmpstr1 = strmid(strsplit(Header(20),/extract),0)
   WaveLengths(ifile) = tmpstr1(1)

   ; Get observation date and time
   tmpstr2 = strmid(strsplit(Header(6),/extract),0)
   ObsDates(ifile)=strmid(tmpstr2(2),1,10)
   tmpstr3 = strmid(strsplit(Header(8),/extract),0)
   ObsFullDates(ifile) = tmpstr3(1)
   
   print, 'EIT wavelength: ', WaveLengths(ifile)
   print, 'Observation time: ',ObsFullDates(ifile)
   
endfor

;check all input images have the same size
;check all wavelengths are different
for ifile=0,nfile-1 do begin
   for jfile=ifile+1,nfile-1 do begin
      if npix(ifile) ne npix(jfile) then begin
         print,'WARNING!! Not all images have the same size. Cannot proceed.'
         retall
      endif
      if WaveLengths(ifile) eq WaveLengths(jfile) then begin
         print,'WARNING!! Two or more input file(s) have the same wavelength!'
         print,'Aborting....'
         retall
      endif
   endfor
endfor
  
; determine if image should be resized
case Nx of
   1024: begin
      if noresize eq 0 then begin
         NxRsz = Nx/2
         NpixRsz = npix(0)/4
         nI = NxRsz
         nJ = NxRsz
         FileDat = 'los_euv_'+ObsDates(0)+'.dat'
         print,''
         print, 'Resizing images to ',nI,' x ',nJ
         ; Assign arrays for output
         DataRsz = fltarr(NpixRsz)
         DataRszAll = fltarr(NpixRsz,nfile)
      endif
   end
   512: begin
      noresize=1
   end
   else: begin
      print, 'WARNING! Image should be either 1024x1024 or 512x512.'
      print, 'Aborting...'
      retall
   end
endcase
if noresize eq 1 then begin 
   nI = Nx
   nJ = Ny
   print,''
   print, 'Keeping original size ',nI,' x ', nJ
   FileDat = 'los_euv_orig_'+ObsDates(0)+'.dat'
   DataAll = fltarr(npix(0),nfile)
endif

; read image data from input fits files 
for ifile=0,nfile-1 do begin
   Data = readfits(filenames(ifile),/silent)
   if noresize eq 0 then begin   
      for i=0L,NxRsz-1 do begin
         for j=0L,NxRsz-1 do begin
            DataRsz(i*NxRsz+j)=(Data(2*(i*Nx+j))+Data(2*(i*Nx+j)+1)+$
                                Data((2*i+1)*Nx+2*j)+Data((2*i+1)*Nx+2*j+1))/4
         endfor
      endfor
     
      ; sort RESIZED data according to increasing wavelength
      if WaveLengths(ifile) eq 171 then DataRszAll(*,0) = DataRsz   
      if WaveLengths(ifile) eq 195 then DataRszAll(*,1) = DataRsz
      if WaveLengths(ifile) eq 284 then DataRszAll(*,2) = DataRsz
      
   endif
   if noresize eq 1 then begin
       ; sort ORIGINAL size data according to increasing wavelength
      if WaveLengths(ifile) eq 171 then DataAll(*,0) = Data
      if WaveLengths(ifile) eq 195 then DataAll(*,1) = Data
      if WaveLengths(ifile) eq 284 then DataAll(*,2) = Data
   endif
endfor

; write output file
print,''
print,'Writing ASCII (TEC)  file ',FileDat
print,''
OutputLine = fltarr(nfile)

openw,lun,FileDat,/get_lun
printf,lun,' TITLE=" EUV Data converted from FITS file "'
printf,lun,'VARIABLES = "X", "Y", "euv171", "euv195", "euv284"'
printf,lun,'ZONE T="LOS Image", I= ',nI,' J= ',nJ,', K=1, F=POINT'

x=-19.8
d=-2*x/(nI-1)
for i=0L,nI-1 do begin
   x= x+d
   y = -19.8
   for j=0L,nI-1 do begin
       y=y+d
       if noresize eq 1 then OutputLine(0:2) = DataAll(i*nI+j,0:2)
       if noresize eq 0 then OutputLine(0:2) = DataRszAll(i*nI+j,0:2)
       printf,lun, format = '(1e14.6,1e14.6,1e14.4,1e14.4,1e14.4)',$
              x,y,OutputLine
    endfor
endfor

print,'Conversion done'

end
