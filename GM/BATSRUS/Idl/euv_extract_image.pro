;------------------------------------------------------------------------------
; euv_extract_image.pro
;------------------------------------------------------------------------------
; Written by C. Downs (2009)
;
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

function euv_extract_image, ImageIn_II, ProjectedX_I, ProjectedY_I, $
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
