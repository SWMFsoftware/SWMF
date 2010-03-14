;==========================================================================
; Written by Rona Oran , 14 Mar 2010
;
; Purpose:
;
;  1. Read EIT SXR fits file
; 
;  2. Resize image to 512 x 512 if needed (same as synthesized images) 
;  
;  3. Write data into a .dat ASCII file, in the same
;  format that synthesized LOS SXR are outputted by BATSRUS (LC component).
;
; Usage:
;   .r sxr_fits_to_ascii
;
;   you will be prompted to enter a fits file name. Files should
;   contain images either of size 1024x1024 or 512x512. If not,
;   routine will be aborted with a warning.
;   
;   noresize=1
;   .r sxr_fits_to_ascii
;
;   This option allows to output the combined SXR data without
;   resizing the image.  This is mostly intended for testing purposes.
;
; Output:
;
;   A file called los_sxr_xxxx.dat (in TECplot format) will be
;   created ('xxxx' stands for the observation date as read from fits file)
;   
;=========================================================================
print, ''
print, 'This routine creates a TEC(.dat) file from SXR fits file.'
print, 'Use euv_fits_to_ascii for EUV images.'
print, ''

; prompt user for filename
print, 'Enter input filename '
askstr,'filename ',filename,1

; read fits file, extract some header info 
Data = readfits(filename, Header,/silent)
print, ''
print, 'Processing fits file ', filename
  
; Get image dimensions
s=size(Data)
Nx=s(1)
Ny=s(2)
npix = Nx*Ny
if Nx ne Ny then begin
   print,'WARNING!! Fits image must be square.'
   print, 'Aborting...'
   retall
endif
; Get EUV wavelength
tmpstr1 = strmid(strsplit(Header(20),/extract),0)
WaveLength = tmpstr1(1)

; Get observation date and time
tmpstr2 = strmid(strsplit(Header(6),/extract),0)
ObsDate = strmid(tmpstr2(2),1,10)
tmpstr3 = strmid(strsplit(Header(8),/extract),0)
ObsFullDate = tmpstr3(1)
   
print, 'EIT wavelength: ', WaveLength
print, 'Observation time: ',ObsFullDate
   
; determine if image should be resized
case Nx of
   1024: begin
      if noresize eq 0 then begin
         NxRsz = Nx/2
         NpixRsz = npix/4
         nI = NxRsz
         nJ = NxRsz
         FileDat = 'los_sxr_'+ObsDate+'.dat'
         print,''
         print, 'Resizing images to ',nI,' x ',nJ
         ; Assign arrays for output
         DataRsz = fltarr(NpixRsz)
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
   FileDat = 'los_sxr_orig_'+ObsDate+'.dat'
endif

if noresize eq 0 then begin   
   for i=0L,NxRsz-1 do begin
      for j=0L,NxRsz-1 do begin
         DataRsz(i*NxRsz+j)=(Data(2*(i*Nx+j))+Data(2*(i*Nx+j)+1)+$
                             Data((2*i+1)*Nx+2*j)+Data((2*i+1)*Nx+2*j+1))/4
      endfor
   endfor
endif

; write output file
print,''
print,'Writing ASCII (TEC)  file ',FileDat
print,''

openw,lun,FileDat,/get_lun
printf,lun,' TITLE=" EUV Data converted from FITS file "'
printf,lun,'VARIABLES = "X", "Y", "sxr"'
printf,lun,'ZONE T="LOS Image", I= ',nI,' J= ',nJ,', K=1, F=POINT'

x=-19.8
d=-2*x/(nI-1)
for i=0L,nI-1 do begin
   x= x+d
   y = -19.8
   for j=0L,nI-1 do begin
       y=y+d
       if noresize eq 1 then Output = Data(i*nI+j)
       if noresize eq 0 then Output = DataRsz(i*nI+j)
       printf,lun, format = '(1e14.6,1e14.6,1e14.4)',x,y,Output
    endfor
endfor

print,'Conversion done'

end
