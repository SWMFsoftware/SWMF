;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;------------------------------------------------------------------------------
;---- collection of routines to download, reduce, and save SOHO EIT images ----
;----  for comparision with los_euv output from the SWMF ----------------------
;
; **** NEED TO COMPILE THIS FILE TWICE IN ORDER FOR THE ROUTINES TO WORK
; **** (the IDL interpreter needs two passes to recognize the custom functions)
; **** (so ignore the error the first time                                    )
;
; routines in this file:
;     get_eit_set       (Main routine, call this from prompt)
;     center_image      (Centers a CCD image)
;     extract_image     (Extracts a FOV, and resizes)
;     hae_to_hgi        (converts HAE coordinates to HGI)
;     rot_matrix        (rotation matrix generator)
;     euler_rot         (does rotation by euler angles)
;     hgi_to_input      (prints PARAM.in input, option to account for SimTime)
;     mk_tv_image       (makes a bytscaled image for display)
;     save_eit_pngs     (saves color pngs of eit images, synth or observations)
;     display_eit_set   (displays color eit images, synth or observations)
;     read_los_tec      (reads an SWMF ascii tecplot file saved by los euv tec)
;------------------------------------------------------------------------------

;==============================================================================
; get_eit_set
;------------------------------------------------------------------------------
; Automated function to use the IDL interface to the virtual solar observatory 
;  (VSO) to obtain, calibrate, and save the nearest set of SOHO EIT 171, 195, 
;  and 284 images to the time supplied.
;
; the Solar SoftWare Idl package (SSW) must be installed, up to date and
;  properly loaded/linked in IDL for this to work properly. 
;  The additional SSW packages needed are SOHO/EIT, and VSO.
;  ALSO, the SSWDB calibration database for SOHO/EIT must be installed and up
;  to date at least for the time requested or the images will NOT calibrate
;  correctly. (eit_prep might also crash if SSWDB is not configured right)
;
; This routine will do the following steps:
;  a) search the VSO using vso_search and find the nearest 3 filter set to a 
;     given input time. 
;
;  b) download the 3 raw UNCALIBRATED Level 0 EIT images [171, 195, 284]
;     using the vso_get interface 
;
;  c) call eit_prep to calibrate the images to units of DN/s.
;
;  d) save the calibrated images to a tecplot ascii file (option, default ON)
;
;  e) color scale the images and display the images (option, default ON)
;
;  f) save color scaled png files (option, default ON)
;
;  g) optionally return the image arrays/Headers/Soho position to IDL variables
;      (need to specify the keywords: Images, Headers, HgiXyz)
;
; **** the default is to downscale the images to match the 512x512 los_euv
;       output of the SWMF. 
;       To prevent significant resolution loss use nPix=1024, rFov=1.4
;
; **** NOTE EIT sets are typically taken every 6 hours at 01, 07, 13, and 19 UT
;       This routine is only for obtaining 3 filter sets. Only the 195 band 
;       takes images at a higher cadence, thus a manual search is probably better 
;       for any time dependent study.
;
; **** Also, this uses the VSO IDL interface in SSW. As of writing this would
;      crash using my 64bit IDL on a linux machine. Using 32 bit mode (idl -32)
;      seems to work fine.
;      
; INPUT: 
;     InDate: a string with the requested time format = 'yyyy/mm/ddThh:mm:ss'
;            (leave the capital T as is), other formats might be unpredictable
;
; OUTPUTS: (files saved to FitsDir, .dat and .pngs to OutDir, 
;           also see keywords for optional output)
;
; KEYWORDS:
;     FitsDir: string with path to directory to save the fits files 
;              (default = ./)
;
;     OutDir: string with path to directory to save the .dat and png files 
;              (default = ./) 
;
;     UseFullPath: use this to make the full eit database subdirectories on top
;                  of the FitsDir path and place the files in the same relative
;                  locations as on the server (FitsDir/eit/lz/YYYY/MM/filename)
;
;     rFov: Field of view radius in solar radii units to save 
;            (default value to match default SWMF saves is 1.98 Rs,
;              note: approx true fov of EIT instrument is ~1.4 Rs)
;
;     nPix: desired number of pixels along each axis for output image 
;            (default value to match default SWMF saves is 512)
;
;     Images: this keyword will recieve a [nPix, nPix, 3] array with the
;             reduced, centered, and resized images in order: 171,195,284
;
;     Headers: this keyword will return the final headers as well 
;               as a [nLine, 3] string array
;
;     ProjectedX_I: return a vector of length values in solar radii units, 
;                    that correspond to the projected position of each pixel 
;                    on the plane interestecting the sun and normal to the 
;                    observer. (assuming square, so same values for x and y)
;
;     HgiXyz: position of the craft in Heliographic Inertial Coordinates
;              This is the coordinate system needed for the los saves in SWMF
;
;     SimDate: Set this to a time string specifying the input simulation time of
;               the SWMF run. If the observations are at a much different time
;               than the one simulated, will need to account for solar rotation
;               This will rotate the HgiXyz position so that the correct 
;               viewing angle is obtained. (Applicable for SWMF runs in HGR)
;               (the same time string format as InDate: 'yyyy/mm/ddThh:mm:ss')
;
;     NoTecFile: Do not save a tecplot ascii file
;
;     NoDisplay: Do not Display the images to the screen
;
;     NoPngs: Do not save the png files
;
;     NoCal: Do Not calibrate the images (FOR TESTING ROUTINE ONLY!)
;
; EXAMPLES:
;   -Default behaivior: 
;     get_eit_set, '1996/08/27T01:00:00'
;
;   -Save the level zero fits to a specific directory head, and the data/pngs
;     to a separate location. Also use full database path for raw data:
;        get_eit_set, '2002/08/24T01:00:00', FitsDir='./soho_data/', $
;                                  OutDir='./output/', /UseFullPath
;
;   -Save the EIT images at full resolution (small shift/stretch is applied)
;        get_eit_set, '2003/04/17T08:00:00', nPix=1024, rFov=1.4
;
;   -Account for the simulation time being different than the observations
;     (will only modify the output Hgi positions)
;        get_eit_set, '2005/01/01T19:00:00', SimDate='2005/01/01T00:00:00'
;
;   -Do not display the images:
;        get_eit_set, '2003/08/27T13:00:00', /NoDisplay
;              
;   -Forego tecplot files/pngs, return the reduced images to Idl keywords:
;        get_eit_set, '1997/12/31T19:00:00', Images=Images, Headers=Headers, $
;                       HgiXyz=HgiXyz, /NoDisplay, /NoTecFile, /NoPngs
;
;   -Forego saving TecFile and PNGs, return the reduced images to Idl keywds:
;     also get the header and the hgi coordinates
;        get_eit_set, '1999/12/31T19:00:00', Images=Images, Headers=Headers, $
;                                  HgiXyz=HgiXyz, /NoTecFile, /NoPngs
;
;   -Test functionality without calling eit_prep for calibration
;     (images won't look right because of lack of calibration)
;        get_eit_set, '2008/03/25T01:00:00', /NoCal
;
;------------------------------------------------------------------------------

pro get_eit_set, InDate, $
   FitsDir = FitsDir, $
   OutDir = OutDir, $
   UseFullPath = UseFullPath, $
   rFov = rFov, $
   nPix = nPix, $
   Images = Images, $
   Headers = Headers, $
   ProjectedX_I = ProjectedX_I, $
   HgiXyz = HgiXyz, $
   SimDate = SimDate, $
   NoTecFile = NoTecFile, $
   NoDisplay = NoDisplay, $
   NoPngs = NoPngs, $
   NoCal = NoCal

;----------------------------------------------------------
;-------- FIRST FIND AND GET THE RAW DATA -----------------
;----------------------------------------------------------
; check requested directories
if not keyword_set(FitsDir) then FitsDir = './'
;check to make sure there is a / on the end (for compatibility w/ later calls)
if strmid(FitsDir,strlen(FitsDir)-1,1) ne '/' then FitsDir=FitsDir+'/'

if not keyword_set(OutDir) then OutDir = './'
;check to make sure there is a / on the end (for compatibility w/ later calls)
if strmid(OutDir,strlen(OutDir)-1,1) ne '/' then OutDir=OutDir+'/'

; use TAI time for time sorting
; (This is a SSW dependency, could make it standalone, but need SSW for VSO anyway...)
InTime = utc2tai(anytim2utc(InDate))

; generate a +- 8hour window from the requested date
TimeStringMinus = strjoin( strsplit( utc2str( tai2utc(InTime - 8*3600.) ), '-', $
                     /extract) , '/')
TimeStringPlus = strjoin( strsplit( utc2str( tai2utc(InTime + 8*3600.) ), '-', $
                     /extract) , '/')

TimeRange = TimeStringMinus + '-' + TimeStringPlus

; get the catolog listing structure for this range
List = vso_search(date=TimeRange,inst='eit')

; get the image times in TAI format
ObsTime_I = utc2tai(List.Time.Start)

; define some variables
nBand = 3
Bands = ['171', '195', '284']
Index_I = intarr(3)
Filenames_I = strarr(3)

; loop over bands, find nearest image set. (find closest 171 then anchor off this)
for iBand=0, nBand-1 do begin

   ; get times for images that match this band and are over 512x512 in size
   ; (a 512x512 image will return a vesorecord.size as 523.0)
   Match = where(List.Wave.Min eq Bands[iBand] and List.Size gt 600., Count)
   ; if found no matches then relax size limit
   if Count lt 1 then Match = where(List.Wave.Min eq Bands[iBand])
   ; if still found no matches then die
   if Count lt 1 then begin
      print,'ERROR!'
      print,'get_eit_set: could not find an image for band: ', Bands[iBand]
      print,'For time range: ', TimeRange
      print,'EXITING'
   endif

   ; these are now the times that match this band
   ObsTimeBand_I = ObsTime_I[Match]

   ; find the closest time to the requested time
   Index_I[iBand] = where(abs(ObsTime_I - InTime) eq min(abs(ObsTimeBand_I - InTime)))

   ; now set the intime to selected 171 image, this way 195 and 284 will be closest
   if iBand eq 0 then InTime = ObsTime_I[Index_I[iBand]]

   ; extract the eit filename from the fileid tag
   ; *** note the archive full path seems to change with date, but last parts
   ;     of the path are constant (eit/lz/YYYY/MM/filename part)
   StringLine = strsplit(List[Index_I[iBand]].fileid,'/', /extract)
   nStrings = n_elements(StringLine)
   FileNames_I[iBand] = StringLine[nStrings-1]

endfor

; this is the structure array of the good images
GoodList = List[Index_I]

print,'Downloading EIT images for times: ', GoodList.Time.Start

; now download them with VSO get (will not re-dowload if complete file exists)
if not keyword_set(UseFullPath) then begin
   mk_dir,FitsDir,/verbose       ; make directory if needed
   ; if not specifying subdirectories, can DL as a vector
   Status = vso_get(GoodList,out_dir=FitsDir)
   FileNames_I = FitsDir + FileNames_I
endif else begin
   ;if using full path, get file location for each and make subdir
   ; this will look like FitsDir/eit/lz/YYYY/MM/filename
   for iBand=0,nBand-1 do begin
      StringLine = strsplit(GoodList[iBand].fileid,'/', /extract)
      nStrings = n_elements(StringLine)
      FileDir = FitsDir + strjoin(StringLine[nStrings-5: nStrings-2],'/') + '/'
      mk_dir,FileDir,/verbose       ; make directory if needed
      Status = vso_get(GoodList[iBand],out_dir=FileDir)
      FileNames_I[iBand] = FileDir + FileNames_I[iBand]
   endfor
endelse

print,''
print,'----- DONE DOWNLOADING IMAGES -----'
print,'EIT 171: ', FileNames_I[0]
print,'EIT 195: ', FileNames_I[1]
print,'EIT 284: ', FileNames_I[2]
print,''

;----------------------------------------------------------
;-------- NOW CALIBRATE THE DATA -----------------------------
;----------------------------------------------------------
if not keyword_set(NoCal) then begin
   eit_prep, FileNames_I, Header_II, Image_III

endif else begin
   ; if not calibrating just read the images
   Image_III = fltarr(1024,1024,nBand)
   Header_II = strarr(77,nBand)
   for iBand=0,nBand-1 do begin
      Image_III[*,*,iBand] = readfits(FileNames_I[iBand], Header)
      Header_II[*,iBand] = Header
   endfor
endelse

;----------------------------------------------------------
;-------- NOW CENTER AND RESIZE ---------------------------
;----------------------------------------------------------
; *** these are defaults to mimic default los_euv saves of SWMF
;     can change these keywords for higher resolution/partial frames
if not keyword_set(rFov) then rFov = 1.98
if not keyword_set(nPix) then nPix = 512

; set this to the square range you want to extract about sun-center
RangeWant = [-rFov, rFov]
; set these to desired pixel resolution
nI = nPix
nJ = nPix
; this is where the aligned and resized images go
ImageNew_III = fltarr(nI, nJ, nBand)

for iBand=0, nBand-1 do begin

   ; center of image in pixel indeces   
   CRPIX1 = float((strsplit(Header_II[37,iBand],/extract))[2])
   CRPIX2 = float((strsplit(Header_II[38,iBand],/extract))[2])
   ; plate scale in arcseconds
   CDELT1 = float((strsplit(Header_II[41,iBand],/extract))[2])
   ; coordinates of craft in Heliocentric Aries Ecliptic (HAE) in km
   HEC_X = float((strsplit(Header_II[49,iBand],/extract))[2])
   HEC_Y = float((strsplit(Header_II[50,iBand],/extract))[2])
   HEC_Z = float((strsplit(Header_II[51,iBand],/extract))[2])
   ; Hae Coords and distance in solar radii units
   Hae_D = [HEC_X, HEC_Y, HEC_Z]/6.96d5
   Distance = sqrt(total(Hae_D^2))

   ; center the image
   ImageCenter_II = center_image(Image_III[*,*,iBand], $
                        [CRPIX1, CRPIX2], CDELT1, Distance, $
                        ProjectedX_I)

   ; interpolate the image to the desired uniform grid in solar radii units
   ImageNew_III[*,*,iBand] = $
      extract_image( ImageCenter_II, ProjectedX_I, ProjectedX_I, $
         nPixOutX=nI, nPixOutY=nJ, $
         xRange=RangeWant, yRange=RangeWant, $
         ProjectedXOut_I=ProjectedXOut_I, ProjectedYOut_I=ProjectedYOut_I)

   ; get and print HGI location information for the first image
   ; (will need to use this as input coords for the SWMF)
   if iBand eq 0 then begin
      Hgi_D = hae_to_hgi(Hae_D)

      ; if SimDate is set, rotate HGI coords so that viewing location
      ; will be the correct line of sight in the simulation
      if keyword_set(SimDate) then begin
         SimTimeSec = utc2tai(anytim2utc(SimDate))
      endif

      ; this routine prints the coords and returns final value to HgiXyz
      hgi_to_input, Hgi_D, ObsTimeSec=InTime, SimTimeSec=SimTimeSec, $
         HgiOut_D = HgiXyz

   endif

endfor

; put the Images and headers into their keywords(will be returned if specified)
Images = ImageNew_III
Headers = Header_II
ProjectedX_I = ProjectedXOut_I

;----------------------------------------------------------
;-------- NOW SAVE A TECPLOT FORMAT DATAFILE --------------
;----------------------------------------------------------

if not keyword_set(NoTecFile) then begin
   ; get the image date time of the first image
   ObsTimeFull = GoodList[0].Time.Start

   ; make the output directory if needed
   mk_dir, OutDir, /verbose
   ; this is the save file name
   FileDat = OutDir + 'los_euv_'+ObsTimeFull+'.dat'
 
   print,''
   print,'Writing ASCII (TEC)  file ',FileDat
   print,''

   openw,lun,FileDat,/get_lun
   printf,lun,' TITLE=" EUV Data converted from FITS file "'
   printf,lun,'VARIABLES = "X", "Y", "euv171", "euv195", "euv284"'
   printf,lun,'ZONE T="LOS Image", I= ',nI,' J= ',nJ,', K=1, F=POINT'


   for i=0L,nI-1 do begin
      for j=0L,nJ-1 do begin
          printf,lun, format = '(1e14.6,1e14.6,1e14.4,1e14.4,1e14.4)',$
                 ProjectedXOut_I[i], ProjectedYOut_I[j], ImageNew_III[i,j,*]
       endfor
   endfor

   free_lun, lun

   print,'Conversion done'

endif

;----------------------------------------------------------
;-------- NOW DISPLAY SOME COLOR PNGs ---------------------
;----------------------------------------------------------
; ** hardcode the png and display images to have a field of view no greater
;    than 1.4. if rFov is larger, reduce png pixel size
Range = min([1.4,rFov])
nPixPng = nPix
if rFov gt 1.4 then nPixPng = nPix*1.4/rFov

; these are the default desired plotranges (in DN/s values)
; Modifying these values here will work for both display and png saving
PlotRange_II= fltarr(2,nBand)
;171
PlotRange_II[0,0] = 0.2
PlotRange_II[1,0] = 3e3
;195
PlotRange_II[0,1] = 0.2
PlotRange_II[1,1] = 3e3
;284
PlotRange_II[0,2] = 0.1
PlotRange_II[1,2] = 1e3

if not keyword_set(NoDisplay) then begin

   ; can change this, but note will be 3x this value in length with all 3
   nPixDisp = 300

   ; call display routine
   display_eit_set, ImageNew_III, ProjectedXOut_I, nPixDisp, $
      rFov=Range, PlotRange_II = PlotRange_II

endif

;----------------------------------------------------------
;-------- NOW SAVE SOME COLOR PNGs ------------------------
;----------------------------------------------------------
; call the image save routine
if not keyword_set(NoPngs) then begin
   
   ; make the output directory if needed
   mk_dir, OutDir, /verbose
   
   ; get the png filenames, its long but will be unique for a given time
   PngFiles_I = OutDir + 'los_euv_'+GoodList[*].Time.Start+'_'+Bands+'.png'

   ; call save routine
   save_eit_pngs, ImageNew_III, ProjectedXOut_I, nPixPng, PngFiles_I, $
      rFov=Range, PlotRange_II = PlotRange_II
endif

print,''
print,'DONE!'
print,''

end
;==============================================================================
; center_image.pro
;------------------------------------------------------------------------------
; Routine to take in an image and center it based on pixel center.
;  Missing edges are filled in with values of 1E-06
;
;  NOTE that this assumes that the image is square (i.e. a normal CCD image)
;   any sub-frame extraction / shrinking if needed is applied LATER
;
;  Also uses angular size seen by pixel, and radial distance
;  to determine physical values of x,y in plane normal to observation vector
;  that intersects the sun center. (y projected along solar north)
;
;  INPUTS:
;     Image_II : 2x2 image array
;     SunCenter_Xy : Pixel location of Sun Center (x,y pair)
;     cDelt : angular size viewed by a pixel (arcseconds)
;     Distance : Distance of instrument from sun (Rsun units)
;
;  OUTPUTS:
;     ImageOut_II : centered image
;     ProjectedX_I: array of x-values corresponding to every pixel column
;
;------------------------------------------------------------------------------

function center_image,Image_II, SunCenter_Xy, cDelt, Distance, $
            ProjectedX_I

; find the number of pixels
nPix = (size(Image_II))[1]

; calculate offsets 
OffsetX = nPix/2 - SunCenter_Xy[0]
OffsetY = nPix/2 - SunCenter_Xy[1]

; setup interpolation arrays (these are pixel indexes)
IndexX_I = findgen(nPix) - OffsetX
IndexY_I = findgen(nPix) - OffsetY

; interpolate the new image to these indeces
ImageOut_II = interpolate(Image_II, IndexX_I, IndexY_I, /grid, missing=1e-6)

; Now want to return the projected locations of the pixels on the plane
; intersecting the sun and normal to observation vector (Rsun units)
; (assuming its a square image --> y positions equal x positions)

; first calculate arc second values
ArcSecX_I = (dindgen(nPix)-nPix/2+0.5)*cDelt

; now go from arcseconds to physical distance on this plane
ProjectedX_I = Distance*sin(ArcSecX_I/3600.*!dtor)

return, ImageOut_II

end

;==============================================================================
; extract_image.pro
;------------------------------------------------------------------------------
; Routine to take in an Image with size and position information 
;  use this to extract and/or resize all or part of the image
;
;  Main use is to reduce the data set size to something manageable, but also
;   helpful for extracting sub-frames.
;
;
;  INPUTS:
;     ImageIn_II : an image array
;     ProjectedX_I : the projected x location of each pixel (row) in Rs units
;     ProjectedY_I : the projected y location of each pixel (column) in Rs
;
;  OUTPUTS:
;     ImageOut_II : the new output image
;
;  KEYWORDS:
;     (inputs)
;        xRange : the x range you want to extract
;        yRange : the y range you want to extract
;        nPixOutX : the desired output image size in the x-direction
;        nPixOutY : the desired output image size in the y-direction
;     (outputs)
;        ProjectedXOut_I : the new values of ProjectedX_I
;        ProjectedYOut_I : the new values of ProjectedY_I
;
;------------------------------------------------------------------------------

function extract_image, ImageIn_II, ProjectedX_I, ProjectedY_I, $
   xRange=xRange, $
   yRange=yRange, $
   nPixOutX=nPixOutX, $
   nPixOutY=nPixOutY, $
   ProjectedXOut_I = ProjectedXOut_I, $
   ProjectedYOut_I = ProjectedYOut_I

SubPrefix = 'extract_image.pro: '

; find the number of pixels of the original image
nPixX = (size(ImageIn_II))[1]
nPixY = (size(ImageIn_II))[2]

; Note, this routine converts to a Regular grid in X,Y. And original images are
; on a uniform grid of Angle ---> ever so slight perspective switch between
; them., --> will need to interpolate the image even if just reducing size.

if not keyword_set(xRange) then xRange = [min(ProjectedX_I), $
      max(ProjectedX_I)]

if not keyword_set(yRange) then yRange = [min(ProjectedY_I), $
      max(ProjectedY_I)]

; if did not specify a pixel size then use ranges to determine size
; this way you do not oversample beyond observation resolution
if not keyword_set(nPixOutX) then begin
   nPixOutX = n_elements(where(ProjectedX_I lt xRange[1] and $
                               ProjectedX_I gt xRange[0], Count))
   if Count lt 2 then begin
      print,SubPrefix+'desired range is TOO small for this image!'+$
            'got nPixOutX = ',nPixOutX
      RETURN, 0
   endif
endif

if not keyword_set(nPixOutY) then begin
   nPixOutY = n_elements(where(ProjectedY_I lt yRange[1] and $
                               ProjectedY_I gt yRange[0], Count))
   if Count lt 2 then begin
      print,SubPrefix+'desired range is TOO small for this image!'+$
            'got nPixOutY = ',nPixOutY
      RETURN, 0 
   endif
endif

; positions of desired REGULAR X,Y grid
PositionsX_I = dindgen(nPixOutX)/(nPixOutX-1)*(xRange[1]-xRange[0]) + xRange[0]
PositionsY_I = dindgen(nPixOutY)/(nPixOutY-1)*(yRange[1]-yRange[0]) + yRange[0]

; interpolate desired positions to pixel coordinate vectors
; these range from 0-nPixX(Y)-1
PixX_I = interpol( dindgen(nPixX), ProjectedX_I, PositionsX_I)
PixY_I = interpol( dindgen(nPixY), ProjectedY_I, PositionsY_I)

; interpolate the image to this grid
ImageOut_II = interpolate(ImageIn_II, PixX_I, PixY_I, /grid, missing=1e-6)

; put in the new values
ImageOut_II = ImageOut_II
ProjectedXOut_I = PositionsX_I
ProjectedYOut_I = PositionsY_I

return, ImageOut_II

end

;==============================================================================
; hae_to_hgi.pro
;------------------------------------------------------------------------------
; Quick routine to convert Hae coordinates to Hgi coordinates
;
;  INPUTS:
;     Hae_D : position in Heliocentric Aries Ecliptic coords (HAE)
;              in Solar Radii units
;
;  OUTPUTS:
;     Hgi_D : position in HelioGraphic Inertial coordinates (HGI)
;              in Solar Radii units
;
;------------------------------------------------------------------------------

function hae_to_hgi,Hae_D

; HAE angles relative to HGI
HaeLongitude = 75.77
HaeInclination = 7.25

Hgi_D = euler_rot(Hae_D, HaeLongitude, HaeInclination, 0.0)

return, Hgi_D

end

;==============================================================================
; rot_matrix.pro
;------------------------------------------------------------------------------
; simple function to return a rotation matrix for rotation about any of the 
;  principle axes (x, y, or z) by an input angle (in radians)
;
; this is entirely based on the rot_matrix fortran routine in 
;  GM/BATSRUS/srcPostProc/CON_geopack_internal.f90 but this time it returns
;  RIGHT handed rotation matrices. 
;    (when rotation definition is: v_D' = matmul(rot_matrix(), v_D) )
;
; INPUT:
;     Dir: integer for which axis to rotatate about (x_=0, y_=1, z_=2)
;     Angle: rotation anlge in radians (right handed)
;
; OUTPUT:
;     returns a 3x3 rotation matrix
;
;------------------------------------------------------------------------------

function rot_matrix,Dir,Angle

;=== main 3x3x3 array that Direction is pulled out of
tt = dblarr(3,3,3)
;x
tt[*,0,0]=[0,0,0]
tt[*,1,0]=[0,0,-1]
tt[*,2,0]=[0,1,0]
;y
tt[*,0,1]=[0,0,1]
tt[*,1,1]=[0,0,0]
tt[*,2,1]=[-1,0,0]
;z
tt[*,0,2]=[0,-1,0]
tt[*,1,2]=[1,0,0]
tt[*,2,2]=[0,0,0]

;== identity matrix
id = dblarr(3,3)
id[*,0]=[1,0,0]
id[*,1]=[0,1,0]
id[*,2]=[0,0,1]

rot_matrix = id      ; make sure its 3x3

if (Angle ne 0.0) then begin
   rot_matrix = cos(Angle)*rot_matrix + sin(angle)*tt[*,*,Dir]
   rot_matrix[Dir,dir] = 1.0
endif

return,rot_matrix

end


;==============================================================================
; euler_rot.pro
;------------------------------------------------------------------------------
; simple function to do rotation about 3 angles given an input vector
;  this is entirely based on the eulerian_matrix fortran routine in 
;  GM/BATSRUS/srcPostProc/CON_geopack_internal.f90
; 
; NOTE! - this does vary slightly from the fortran version but returns the 
;          exact same result as Result = matmul(eulerian_matrix(args), vector)
;
;          the minus sign on the angles is because for some reason
;          in CON_geopack_internal the rot_matrix function returns
;          LEFT handed rotation matrices when you use the form
;          v_D' = matmul(rot_matrix_DD, v_D). But my IDL rot_matrix
;          version is set to return right handed rotations (since it is
;          used as a general function for other things)
;
;  NOTE: fortran matrix multiplication command: V_D' = matmul(M_DD, V_D)
;         is equivalent to IDL command: V_D' = M_DD ## V_D
;         which is the same as you would work it out on paper
;
; INPUTS:
;        Xyz_D: three element vector with x, y, z coordinates
;        InOmega: first Eulerian angle  (in degrees)
;        InTheta: second Eulerian angle (in degrees)
;        InPsi:   third Eulerian angle  (in degrees)
;       
; OUTPUTS:
;        returns a rotated Xyz vector
;------------------------------------------------------------------------------

function euler_rot, Xyz_D, InOmega,InTheta,InPsi

Euler_DD = dblarr(3,3)
x_ = 0
y_ = 1
z_ = 2

;==== convert to radians here!
Omega = InOmega*!dtor
Theta = InTheta*!dtor
Psi = InPsi*!dtor

Euler_DD = rot_matrix(x_, -Theta)##rot_matrix(z_, -Omega)

Euler_DD = rot_matrix(z_, -Psi)##Euler_DD

return, Euler_DD##Xyz_D

end

;------------------------------------------------------------------------------
; hgi_to_input.pro
;------------------------------------------------------------------------------
; Quick routine to take HGI coordinates and print the position input lines 
;  for the los euv part of the #SAVEPLOT call in the BATSRUS PARAM.in file.
;  
; ALSO, because the time of observations one wants to compare to are not 
;  always necessarily the same time as input into the simulation, have optional
;  simulation time and observation time as keyword input. This will add a
;  rotation of the HGI coords so that they now correspond to the correct
;  viewing angle.
;
; *** it is easy to convert a date to TAI seconds!
;     (e.g. SimTimeSec = utc2tai(anytim2utc('2003/08/23T01:00:00'))
;
;  INPUTS:
;     HgiIn_D : position in HelioGraphic Inertial coordinates (HGI)
;              in Solar Radii units
;
;  OUTPUTS:
;     None (prints the output, also can specify HgiOut_D keyword)
;
;  KEYWORDS:
;     HgiOut_D : Optional, return the final HGI coords
;     ObsTimeSec : Time in Tai seconds (double)
;     SimTimeSec : If the simulation time is diferent, rotate the HGI position
;        so that it will correspond the to right viewing position in BATSRUS
;        for the actual observation time.
;
;------------------------------------------------------------------------------

pro hgi_to_input, HgiIn_D, $
   HgiOut_D = HgiOut_D, $
   ObsTimeSec = ObsTimeSec, $
   SimTimeSec = SimTimeSec

x_ = 0
y_ = 1
z_ = 2

Hgi_D = HgiIn_D

RotateString = ''

if keyword_set(SimTimeSec) and keyword_set(ObsTimeSec) then begin

   ; define the same rotation rate variables as in the SWMF
   OmegaCarrington = 2d0*!dPi/(25.38D0*24*3600)

   ; will be positive if SimTime in the future
   TimeDiff = SimTimeSec - ObsTimeSec

   ; just moving HGI here, but its relative to time of the simulation.
   ;  Think of corotating line projected radially outward from the sun.
   ;  If SimTime is in the future then this line will have rotated
   ;  counter clockwise along z-axis (positive rotation) in HGI coords
   ;  over the time difference.
   ; 
   ;  i.e. if SimTime is at say 18:00:00 for some date but the obstime is
   ;  at 03:00:00 for the same day, want the HGI position to be rotated the 
   ;  angular equivalent of fifteen hours of POSITIVE rotation (ahead)
   ;   
   ; Therefore we want a:
   ;  positive angle (counterclockwise rotation) if SimTime is in the future
   ;  negative angle (clockwise rotation) if SimTime is in the past

   Angle = ((OmegaCarrington*(TimeDiff)) mod  (2d0*!dPi))

   Rotate_DD = rot_matrix(z_, Angle)

   Hgi_D = Rotate_DD ## Hgi_D

   RotateString = ' (ROTATED for SimDate: ' +  utc2str(tai2utc(SimTimeSec)) +$
                      ')'

endif

Format = '(F12.2)'

print,''
print, ' HGI locations for SOHO EIT ', RotateString
print,strtrim(string(Hgi_D[x_],format=format),2),'          ObsPosX'
print,strtrim(string(Hgi_D[y_],format=format),2),'          ObsPosY'
print,strtrim(string(Hgi_D[z_],format=format),2),'          ObsPosZ'
print,''

; set output keyword
HgiOut_D = Hgi_D

end

;==============================================================================
; mk_tv_image.pro
;------------------------------------------------------------------------------
; Routine to take in a Image Structure type outputs a BytScled
;  image for display. This also sets the specified color table
;
;  INPUTS:
;     Image_II : a [nPixX. nPixY] image array
;     ProjectedX_I : the projected x location of each pixel (row) in Rs units
;     ProjectedY_I : the projected y location of each pixel (column) in Rs
;     PlotMin : minimum value for scaling
;     PlotMax : maximum value for scaling
;     iColorTable : color table you want to be set
;     nPixOutX : desired number of output pixels in the x-direction
;     nPixOutY : desired number of output pixels in the y-direction
;
;  OUTPUTS:
;     ImageOut_II : Byte Scaled image
;
;  KEYWORDS:
;     xRange : 2 element x Range in solar radii units you want to plot
;     yRange : 2 element y Range in solar radii units you want to plot
;     NoScale : use this if you do not want to return a byte scaled image
;               (but will interpolate it to the correct size/position)
;
;------------------------------------------------------------------------------

function mk_tv_image, Image_II, ProjectedX_I, ProjectedY_I, PlotMin, PlotMax, $
            iColorTable, nPixOutX, nPixOutY, $
            xRange=xRange, $
            yRange=yRange, $
            NoScale=NoScale

; find the number of pixels of the original image
nPixX = (size(Image_II))[1]
nPixY = (size(Image_II))[2]

; Note, this routine converts to a Regular grid in X,Y. And original images are
; on a uniform grid of Angle ---> ever so slight perspective switch between
; them, --> will need to interpolate the image even if plotting the whole
; range.
; NOTE, some processed images will already be on a uniform x,y grid but still
;  need to interpolate if there is any size change, so this is not a large e
;  efficiency loss....

if not keyword_set(xRange) then xRange = [min(ProjectedX_I), $
   max(ProjectedX_I)]

if not keyword_set(yRange) then yRange = [min(ProjectedY_I), $
   max(ProjectedY_I)]

; positions of desired REGULAR X,Y grid
PositionsX_I = dindgen(nPixOutX)/(nPixOutX-1)*(xRange[1]-xRange[0]) + xRange[0]
PositionsY_I = dindgen(nPixOutY)/(nPixOutY-1)*(yRange[1]-yRange[0]) + yRange[0]

; interpolate desired positions to pixel coordinate vectors
; these range from 0-nPixX(Y)-1
PixX_I = interpol( dindgen(nPixX), ProjectedX_I, PositionsX_I)
PixY_I = interpol( dindgen(nPixY), ProjectedY_I, PositionsY_I)

; interpolate the image to this grid
ImageOut_II = interpolate(Image_II, PixX_I, PixY_I, /grid, missing=1e-6)

; if just want the interpolated image, return it
if keyword_set(NoScale) then return, ImageOut_II

; now bytscl the image
ImageOut_II = alog10(ImageOut_II>PlotMin<PlotMax)
ImageOut_II = bytscl(ImageOut_II,min=alog10(PlotMin),max=alog10(PlotMax))

; load the color table (use eit_colors since SSW is needed anyway)
eit_colors,iColorTable, /silent

return, ImageOut_II

end

;==============================================================================
; save_eit_pngs
;------------------------------------------------------------------------------
; Simple routine to take in a set of 3 filter EIT images and output png files
;  Works with synth or observations
;  if the images are observations ,they must be correctly calibrated to DN/s 
;  or they will look very strange
;
; This one ASSUMES that the images are SQUARE and CENTERED for simplicity!!! 
; (though could easily modify call to mk_tv_image for non-square ranges)
;
; INPUTS: 
;     Image_III: a [nPix, nPix, 3] array with eit 171, 195, and 284 images
;     ProjectedX_I: a vector of length values in solar radii units, corresponds 
;                    to the projected position of each pixel on the plane 
;                    interestecting the sun and normal to the observer
;     nPixPng: the pixel size of the output png (eg 512, 1024)
;     PngFiles_I: a three element string array with the desired out filenames
;     
; OUTPUTS: (None to IDL but writes to disk, optionally displays image)
;
; KEYWORDS:
;     rFov: Field of view radius in solar radii units to display
;     PlotRange_II: a [2,nBand] element array with plotting min max values for
;                    each band. Use default values if not specified.
;
;------------------------------------------------------------------------------

pro save_eit_pngs, Image_III, ProjectedX_I, nPixPng, PngFiles_I, $
      rFov=rFov, $
      PlotRange_II=PlotRange_II

nBand=3

; these are the default desired plotranges (in DN/s values)
if not keyword_set(PlotRange_II) then begin
   PlotRange_II= fltarr(2,nBand)
   ;171
   PlotRange_II[0,0] = 0.2
   PlotRange_II[1,0] = 3e3
   ;195
   PlotRange_II[0,1] = 0.2
   PlotRange_II[1,1] = 3e3
   ;284
   PlotRange_II[0,2] = 0.1
   PlotRange_II[1,2] = 1e3
endif

; Color table indexes for eit_colors (blue, green, yellow)
ColorIndexes_I = [42, 43, 44] 

print,''
print,'Now Saving PNG files:'
print,'PNG 171: ', PngFiles_I[0]
print,'PNG 195: ', PngFiles_I[1]
print,'PNG 284: ', PngFiles_I[2]
print,''


; loop over bands, save the images
for iBand=0, nBand-1 do begin

   ; This makes the bytscaled image (yes, the call is horribly long, but it 
   ; can handle more general sizes/distance ranges)   
   ImScl_II = mk_tv_image(reform(Image_III[*,*,iBand]), $
      ProjectedX_I, ProjectedX_I, $
      PlotRange_II[0,iBand], PlotRange_II[1,iBand], $
      ColorIndexes_I[iBand], $
      nPixPng, nPixPng, $
      xRange=[-rFov,rFov], yRange=[-rFov,rFov])

   ; now write the png file
   tvlct,R,G,B,/get
   ImRgb=bytarr(3, nPixPng, nPixPng)
   ImRgb[0,*,*] = R[ImScl_II]
   ImRgb[1,*,*] = G[ImScl_II]
   ImRgb[2,*,*] = B[ImScl_II]
   write_png, PngFiles_I[iBand], ImRgb

endfor

end

;==============================================================================
; display_eit_set
;------------------------------------------------------------------------------
; Simple routine to take in a set of 3 filter EIT images and display them
;  Works with synth or observations
;  
; if the images are observations ,they must be correctly calibrated to DN/s 
;  or they will look very strange
;
; This one ASSUMES that the images are SQUARE and CENTERED for simplicity!!! 
; (though could easily modify call to mk_tv_image for non-square ranges)
;
; INPUTS: 
;     Image_III: a [nPix, nPix, 3] array with eit 171, 195, and 284 images
;     ProjectedX_I: a vector of length values in solar radii units, corresponds 
;                    to the projected position of each pixel on the plane 
;                    interestecting the sun and normal to the observer
;     nPix: the pixel size of one of the output images
;     
; OUTPUTS: (None, but displays images)
;
; KEYWORDS:
;     rFov: Field of view radius in solar radii units to display
;     PlotRange_II: a [2,nBand] element array with plotting min max values for
;                    each band. Use default values if not specified.
;
;------------------------------------------------------------------------------

pro display_eit_set, Image_III, ProjectedX_I, nPix, $
      rFov=rFov, $
      PlotRange_II=PlotRange_II

nBand=3

; these are the default desired plotranges (in DN/s values)
if not keyword_set(PlotRange_II) then begin
   PlotRange_II= fltarr(2,nBand)
   ;171
   PlotRange_II[0,0] = 0.2
   PlotRange_II[1,0] = 3e3
   ;195
   PlotRange_II[0,1] = 0.2
   PlotRange_II[1,1] = 3e3
   ;284
   PlotRange_II[0,2] = 0.1
   PlotRange_II[1,2] = 1e3
endif

; Color table indexes for eit_colors (blue, green, yellow)
ColorIndexes_I = [42, 43, 44] 

; loop over bands, save the images
for iBand=0, nBand-1 do begin

   ; This makes the bytscaled image (yes, the call is horribly long, but it 
   ; can handle more general sizes/distance ranges)   
   ImScl_II = mk_tv_image(reform(Image_III[*,*,iBand]), $
      ProjectedX_I, ProjectedX_I, $
      PlotRange_II[0,iBand], PlotRange_II[1,iBand], $
      ColorIndexes_I[iBand], $
      nPix, nPix, $
      xRange=[-rFov,rFov], yRange=[-rFov,rFov])

   ; first make a window
   ;if (iBand eq 0) and (!d.Name eq 'X') then begin
   if iBand eq 0 then begin
      WinX = nPix*3
      WinY = nPix
      window, /free, xsize=WinX, ysize=WinY
   endif
   ; tv the image, use extra parameters so it should work for ps device too
   tv, congrid(ImScl_II, WinY, WinY), 0.0+iBand/3., 0.0, $
      xsize=1/3., ysize=1.0, /normal

endfor

end

;------------------------------------------------------------------------------
; read_los_tec.pro
;------------------------------------------------------------------------------
; Routine to read in a BATSRUS synthetic los euv image file in tecplot format.
;  The input files need to be saved in the ascii tecplot format.
;
;  the number of pixels and number of images are read from the header
;
; Inputs:
;     FileName : The full file path of a synthetic image file
;
; Outputs:
;     Image_III : output images in format [nPix, nPix, nImages]
;                 (will only be nPix,nPix if only one image)
;
; Keywords:
;     ProjectedX_I : the projected x locations in Rs units in the image plane
;     ProjectedY_I : the projected y locations in Rs units in the image plane
;	
;------------------------------------------------------------------------------

pro read_los_tec, FileName, Image_III,  $
   ProjectedX_I=ProjectedX_I, $
   ProjectedY_I=ProjectedY_I

; set the number of header lines here
nHeader = 3
Header = strarr(nHeader)

; open the file and read in the header
openr, lun, FileName, /get_lun
readf, lun, Header

; count the number of entries on the variables line to determine if it is an euv
; or an sxr file. If header output style changes, this must also change!!!!
nTotal = n_elements(strsplit(header[1],/extract))
nImage = nTotal -4

; extract the third header line to determine pixel size (assuming square)
; If header output style changes, this must also change!!!!
nPix = long((strsplit(header[2],/extract))[4])

; now read in the data, this is a vector of nPix^2 rows and nImage+2 colums
DataVector = fltarr(nImage+2, nPix^2)
readf, lun, DataVector
;close, lun
free_lun, lun

; now put it into useable arrays

Image_III = fltarr(nPix, nPix, nImage)
ImageX = fltarr(nPix,nPix)
ImageY = fltarr(nPix,nPix)

for iPix=0L, nPix-1 do begin
   for jPix=0L,nPix-1 do begin

      ; save x and y values
      ImageX[iPix, jPix] = Datavector[0,nPix*iPix+jPix]
      ImageY[iPix, jPix] = Datavector[1,nPix*iPix+jPix]

      ; save image values
      for iImage=0, nImage-1 do begin
         Image_III[iPix, jPix, iImage] = Datavector[2+iImage,nPix*iPix+jPix]
      endfor

   endfor
endfor

; grids are uniform
ProjectedX_I = ImageX[*,0]
ProjectedY_I = ImageY[0,*]

end
