;--------------------------------------------------------------------------
; euv_fits_to_ascii
;--------------------------------------------------------------------------
;==========================================================================
; Written by Rona Oran , 14 Mar 2010
; Revision history:
; 13 Apr 2010 - 1. Use external functions by Cooper Downs:
;                  euv_center_image.pro :: 
;                  Centers image, fill in edges with 1E-06

;                  euv_extract_image.pro ::
;                  Modifies X-Y pixel positions so that pixels in
;                  fits image (defined on an angular grid) are projected
;                  onto a flat grid lying on a plane though the sun
;                  center ( which is the same as the BATSRUS grid used for
;                  synthesized EUV LOS images).
;         
;               2. Allow user to use default input filenames instead
;                  of typing them in at the prompt.
; Purpose:
;
;  1. Read 3 EIT EUV fits files (corresponding to 3 wavelengths) and merge
;  into one data set
; 
;  2. Resize image to 512 x 512, correct for instrument FOV  to
;  produce an image that can be compared to synthesized LOS images. 
;  
;  3. Write data into a .dat ASCII file, in the same
;  format that synthesized LOS EUV are outputted by BATSRUS (LC component).
;
; Usage:
;
;   .r euv_fits_to_ascii
;
;   INPUT FILES:
;   Make sure that three euv fits files are present in the directory,
;   and that they each correspond to a different wavelength. Otherwise
;   the routine will abort with an error.
;
;   INPUT FILENAMES:
;   You will be asked if you want to use default input file names.
;   This option saves you the need of retyping the files names at
;   every run.
;   In order to use this option, answer "y". Note that you must rename
;   the input files to the default ones prior to using this routine.
;   DEFAULT : fits1,fits, fits2.fits, fits3.fits.
;   It does not matter which file is renamed  1, 2, or 3. The
;   routine will sort files according to wavelength automatically.
;
;   If you answer "n", you can manualy type in the filenames at the prompt.
;
;   KEEPING THE ORIGIAL SIZE:   
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
filename = ''
askstr, 'Use default input filenames(y/n)',filename,1
if filename eq 'y' then begin
   filenames(0) = 'fits1.fits'
   filenames(1) = 'fits2.fits'
   filenames(2) = 'fits3.fits'
   print, ' Make sure three files named fitsX.fits (X = 1,2,3) exist.'
   print,''
endif
if filename eq 'n' then begin
   for ifile=0,nfile-1 do begin
      print, 'Enter input filename ',ifile
      askstr,'filename ',filename,1
      filenames(ifile)=filename
      filename=''
   endfor
endif

; Assign arrays for input
npix=fltarr(nfile)
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
         ;DataRszAll = fltarr(NpixRsz,nfile)
         DataAll = fltarr(NpixRsz,nfile)
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
   Data = readfits(filenames(ifile), Header, /silent)

   ; get header info related to positions
   ; offset with respect to sun center
   CRPIX1 = float((strsplit(Header[37],/extract))[2])
   CRPIX2 = float((strsplit(Header[38],/extract))[2])
   ; angular size of a pixel FOV
   CDELT1 = float((strsplit(Header[41],/extract))[2])
   ; Instrument position
   HEC_X = float((strsplit(Header[49],/extract))[2])
   HEC_Y = float((strsplit(Header[50],/extract))[2])
   HEC_Z = float((strsplit(Header[51],/extract))[2])
   Distance = sqrt(total([HEC_X, HEC_Y, HEC_Z]^2))/6.96d5

   ; center the image
   ImageCenter_II = euv_center_image(Data, [CRPIX1, CRPIX2], $
                        CDELT1, Distance, ProjectedX_I)

   ; set this to the square range you want to extract about sun-center
   RangeWant = [-1.98, 1.98]

   ; interpolate the image to the desired uniform grid 
   ImageNew_II = euv_extract_image( ImageCenter_II, ProjectedX_I, ProjectedX_I, $
      nPixOutX=nI, nPixOutY=nJ, $
      xRange=RangeWant, yRange=RangeWant, $
      ProjectedXOut_I=ProjectedXOut_I, ProjectedYOut_I=ProjectedYOut_I)

   if WaveLengths(ifile) eq 171 then Index = 0
   if WaveLengths(ifile) eq 195 then Index = 1
   if WaveLengths(ifile) eq 284 then Index = 2
  
   ; make sure the pixels are going to the correct indexes in the DataAll array
   for i=0,nI-1 do begin 
      for j=0,nJ-1 do begin
         DataAll[i*nI+j,Index] = ImageNew_II[i,j]
      endfor
   endfor
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

for i=0L,nI-1 do begin
   for j=0L,nJ-1 do begin
       OutputLine(0:2) = DataAll(i*nI+j,0:2)
       printf,lun, format = '(1e14.6,1e14.6,1e14.4,1e14.4,1e14.4)',$
              ProjectedXOut_I[i], ProjectedYOut_I[j], OutputLine
    endfor
endfor

free_lun, lun

print,'Conversion done'

end
