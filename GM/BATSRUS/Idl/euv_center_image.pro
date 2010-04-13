;------------------------------------------------------------------------------
; center_image.pro
;------------------------------------------------------------------------------
; Written by C. Downs (2009)
;
;  Routine to take in an image and center it based on pixel center.
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

function euv_center_image,Image_II, SunCenter_Xy, cDelt, Distance, $
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
