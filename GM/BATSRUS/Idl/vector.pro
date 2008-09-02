;^CFG COPYRIGHT VAC_UM
;=============================================================================
PRO STREAMLINE,U,V,XX,YY,NVECS=nvecs,X0=x0,$
                  XXOLD=xxold,YYOLD=yyold,TRIANGLES=triangles,SEED=seed,$
                  NOERASE=noerase,WHITE=white

;+
; NAME:
;	STREAMLINE
;
; PURPOSE:
;	Draw a vector field with streamlines following the field.
;	Also works with irregular grids using triangulation.
;
; CATEGORY:
;	Graphics, two-dimensional.
;
; CALLING SEQUENCE:
;	STREAMLINE, U, V [, XX, YY]
;
; INPUTS:
;	U:	The X component at each point of the vector field.  This 
;		parameter must be a 2D array.
;
;	V:	The Y component at each point of the vector field.  This 
;		parameter must have the same dimensions as U.
;
; OPTIONAL INPUT PARAMETERS:
;       XX:     X-coordinate array (2D), default is uniform grid with dx=1.
;
;       YY:     Y-coordinate array (2D), default is uniform grid with dy=1.
;
; KEYWORD PARAMETERS:
;	NVECS:	The number of vectors (arrows) to draw.  Default is 200.
;
;	X0:	Positions of the arrows. Default is random.
;
;	SEED: Random number generator seed for proper randomization when VECTOR
;	      is called repaetedly. It has to be a NAMED VARIABLE.
;
;	NOERASE: NOERASE=1 or /NOERASE does not allow erase for VECTOR.
;
;       WHITE:   WHITE=1 or /WHITE forces the vectors to be white.
;
;	The following 3 keyword parameters should be supplied together
;	when the XX,YY parameters describe a nonuniform grid. Do not initialize
;	them! Their purpose is to save on the triangulation time for the next
;	call of VECTOR, if the grid does not change.
;
;	XXOLD:  XX in previous call, on output XX in this call.
;
;	YYOLD:  YY in previous call, on output YY in this call.
;
;       TRIANGLES: Array to store the triangles from the triangulization.
;
; OUTPUTS:
;	A velocity field graph is drawn on the current graphics device.
;
;	If X0 is given it returns the position for the next set of arrows.
;	If XXOLD,YYOLD,TRIANGLES are given, they return the old grid
;	coordinates and the triangles from their triangulation.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	A plot is drawn on the current graphics device.
;
; RESTRICTIONS:
;	Uses procedure vector (see below).
;
; PROCEDURE:
;	Uses vector to draw the arrow heads first, then an very long
;       arrow without a head in the positive direction and then 
;       an arrow without a head in the negative direction.
;       In effect stream lines are drawn with an arrow showing the orientation.
;	If X0 is provided it is used for the location of the arrows instead
;	of random positions.
;
; MODIFICATION HISTORY:
;       05/16/07- G. Toth extracted from existing code.
;-

; normalize vectors
norm = sqrt(u^2+v^2+1.e-30) & u1 = u/norm & v1 = v/norm

; draw arrows
vector,u1,v1,xx,yy,NVECS=nvecs,MAXVAL=1.,$
  NSTEP=6,LENGTH=0.012,HEAD=0.5,$
  XXOLD=xxold,YYOLD=yyold,TRIANGLES=triangles,$
  DYNAMIC=0,SEED=seed,X0=x0,$
  NOERASE=noerase,WHITE=white

; draw streamline along u1;v1
vector,u1,v1,xx,yy,NVECS=nvecs,MAXVAL=1.,$
  NSTEP=1000,LENGTH=2.,HEAD=0.,$
  XXOLD=xxold,YYOLD=yyold,TRIANGLES=triangles,$
  DYNAMIC=0,SEED=seed,X0=x0,$
  /NOERASE,WHITE=white

; draw streamline in the other direction
u1=-u1 & v1=-v1
vector,u1,v1,xx,yy,NVECS=nvecs,MAXVAL=1.,$
  NSTEP=1000,LENGTH=2.,HEAD=0.,$
  XXOLD=xxold,YYOLD=yyold,TRIANGLES=triangles,$
  DYNAMIC=0,SEED=seed,X0=x0,$
  /NOERASE,WHITE=white

end
;=============================================================================
PRO VECTOR,U,V,XX,YY,NVECS=nvecs,MAXVAL=maxval,LENGTH=length,HEAD=head,$
		  NSTEPS=nsteps,X0=x0,DYNAMIC=dynamic,NOERASE=noerase,$
		  XXOLD=xxold,YYOLD=yyold,TRIANGLES=triangles,SEED=seed,$
                  WHITE=white
;+
; NAME:
;	VECTOR
;
; PURPOSE:
;	Draw a velocity (flow) field with arrows following the field 
;	proportional in length to the field strength.  Arrows are composed 
;	of a number of small segments that follow the streamlines.
;	When called repeatedly the arrows can move along the streamlines.
;	Also works with irregular grids using triangulation.
;
; CATEGORY:
;	Graphics, two-dimensional.
;
; CALLING SEQUENCE:
;	VECTOR, U, V [, XX, YY]
;
; INPUTS:
;	U:	The X component at each point of the vector field.  This 
;		parameter must be a 2D array.
;
;	V:	The Y component at each point of the vector field.  This 
;		parameter must have the same dimensions as U.
;
; OPTIONAL INPUT PARAMETERS:
;       XX:     X-coordinate array (2D), default is uniform grid with dx=1.
;
;       YY:     Y-coordinate array (2D), default is uniform grid with dy=1.
;
; KEYWORD PARAMETERS:
;	NVECS:	The number of vectors (arrows) to draw.  Default is 200.
;
;	MAXVAL:	The value of the highest velocity that can possibly occur.
;		Default is calculated as max(sqrt(U^2+V^2)).
;
;	LENGTH:	The length of the arrow corresponding to MAXVAL as a fraction 
;		of the diagonal length of the X-Y box. The default is 0.04.
;
;	HEAD:	The size of the arrow head relative to LENGTH. Default is 0.3.
;
;	NSTEPS:	The number of line_segments forming each arrow. Default is 5.
;
;	X0:	Positions of the arrows. Default is random.
;
;       DYNAMIC:The number of segments by which X0 is shifted. Default is 0.
;
;	SEED: Random number generator seed for proper randomization when VECTOR
;	      is called repaetedly. It has to be a NAMED VARIABLE.
;
;	The following 3 keyword parameters should be supplied together
;	when the XX,YY parameters describe a nonuniform grid. Do not initialize
;	them! Their purpose is to save on the triangulation time for the next
;	call of VECTOR, if the grid does not change.
;
;	XXOLD:  XX in previous call, on output XX in this call.
;
;	YYOLD:  YY in previous call, on output YY in this call.
;
;       TRIANGLES: Array to store the triangles from the triangulization.
;
;	NOERASE: NOERASE=1 or /NOERASE does not allow erase for VECTOR.
;
;       WHITE:   WHITE=1 or /WHITE forces the vectors to be white.
;
; OUTPUTS:
;	A velocity field graph is drawn on the current graphics device.
;
;	If X0 is given it returns the position for the next set of arrows.
;	If XXOLD,YYOLD,TRIANGLES are given, they return the old grid
;	coordinates and the triangles from their triangulation.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	A plot is drawn on the current graphics device.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	NVECS random points within (X,Y) are selected. For each "shot" the 
;	(U,V) field (as bilinearly interpolated) at each point is followed 
;	with NSTEPS line segments. An arrow head is drawn at the end.
;	The length of the arrow corresponding to a velocity MAXVAL is LENGTH.
;	If X0 is provided it is used for the location of the arrows instead
;	of random positions and X0 is advansed by DYNAMIC segments.
;
; MODIFICATION HISTORY:
;	12/2/92	- modified to handle !p.multi (jiy-RSI)
;	          Neal Hurlburt, April, 1988.
;       7/12/94 HJM - Fixed error in weighting factors in function
;                     vel_mybi() which produced incorrect velocity vectors.
;       30/5/96 - R. Keppens changed to allow for dynamic velocity plot.
;		  Added X and Y parameters for correct axis labeling.
;                 X0 stores starting positions, DYNAMIC the index for shifts,
;                 MAXVAL the common normalization for multiple calls.
;       03/6/96 - G. Toth changed a few more things. 
;                 No TITLE keyword parameter, just use !p.title.
;                 The use of !radeg is corrected (used to be the inverse).
;                 Simplified ARRHEAD and ARROW, and merged both into VECTOR.
;		  Aspect ratio is taken into account for drawing arrow heads.
;		  Use NOCLIP=0 for clipping arrows outside box.
;                 Replaced vel_mybi() with IDL's interpolate() function.
;		  LENGTH is made to be relative to the diagonal of the box.
;		  HEAD parameter added for arrowhead size. 0 for streamlines.
;		  NSTEPS now means the number of segments (instead of points).
;	25/11/96- G. Toth added /NOERASE keyword
;	02/12/96- G. Toth and R. Keppens fixed the plotting error for
;		  differrent number of grid points in the X and Y directions.
;       03/27/97- G. Toth and D. Innes generalized VECTOR for non-uniform grids
;	04/07/00- G. Toth added second order integration useful for accurate
;		  stream line drawing. It is used when NSTEPS>5
;       10/01/02- G. Toth added error messages if U and V arrays are wrong
;                 G. Toth added support for irregular grids stored as
;                 a big linear array (ny=1).
;       05/16/07- G. Toth added /WHITE keyword
;       09/02/08- G. Toth use !x.range and !y.range info to limit grid
;-

;Return to caller if an error occurs
;
;on_error,2

if n_elements(U) eq 0 then begin
   print,"Error in vector: first component of the vectorfield is missing"
   retall
endif

if n_elements(V) eq 0 then begin
   print,"Error in vector: second component of the vectorfield is missing"
   retall
endif

if n_elements(U) ne n_elements(V) then begin
   print,"Error in vector: first and second components of the vectorfield"
   print,"                 have different number of elements"
   retall
endif

if max(abs(size(U)-size(V))) ne 0 then begin
   print,"Error in vector: first and second components of the vectorfield"
   print,"                 have different shapes"
   retall
endif

sizes=size(U)
if sizes(0) eq 1 then begin
  ; U is a 1D array (possibly an irregular grid)
  irregular = 1
  nx = sizes(1)
  ny = 1
  if n_elements(XX) eq 0 or n_elements(YY) eq 0 then begin
     print,'Error in vector: one dimensional array without X and Y coordinates'
     retall
  endif
endif else begin
  ; U is a 2D array, assume regular
  irregular = 0
  nx = sizes(1)
  ny = sizes(2)
  ; Default values for optional coordinates
  if n_elements(XX)       eq 0 then xx      = [0,nx-1]
  if n_elements(YY)       eq 0 then yy      = [0,ny-1]
endelse

; Default values for optional parameters
;
if n_elements(NVECS)    eq 0 then nvecs   = 200
if not keyword_set(MAXVAL)   then maxval  = max(sqrt(u^2+v^2))
if n_elements(LENGTH)   eq 0 then length  = 0.04
if n_elements(HEAD)     eq 0 then head    = 0.3
if n_elements(NSTEPS)   eq 0 then nsteps  = 5
if n_elements(DYNAMIC)  eq 0 then dynamic = 0

; Derived parameters and derived defaults
;
xmin=min(XX) & xmax=max(XX)
if !x.range(0) ne !x.range(1) then begin
   xmin = max([xmin, min(!x.range)]) & xmax = min([xmax, max(!x.range)])
endif
ymin=min(YY) & ymax=max(YY)
if !y.range(0) ne !y.range(1) then begin
   ymin = max([ymin, min(!y.range)]) & ymax = min([ymax, max(!y.range)])
endif
xdel=xmax-xmin & ydel=ymax-ymin

if dynamic lt 0      then dynamic=0
if dynamic gt nsteps then dynamic=nsteps

; Check if X0 is provided
random = n_elements(X0) ne 2*nvecs

if random and keyword_set(X0) then print, $
  'VECTOR: Initial position array X0 has incorrect size.',$
  ' Using random positions.'

if random then begin
   ; Random positions within (xmin-xmax,ymin-ymax) if x0 is not defined
   ;
   x0=fltarr(nvecs,2)
   x0(*,0)=xmin+xdel*randomu(seed,nvecs)
   x0(*,1)=ymin+ydel*randomu(seed,nvecs)
endif else begin
   ; Put the x0 arrow positions inside the box by applying a mod function
   ;
   x0(*,0)=x0(*,0)-xdel*floor((x0(*,0)-xmin)/xdel)
   x0(*,1)=x0(*,1)-ydel*floor((x0(*,1)-ymin)/ydel)
endelse

; For a presumably regular grid, check if it is regular or not
if not irregular then $
   ; Check if there is a grid passed in XX, YY
   if n_elements(xx) gt 2 and n_elements(yy) gt 2 then $
      ; Check if the grid is uniform by calculating maximum grid spacing
      if min(xx(1:nx-1,*)-xx(0:nx-2,*)) lt xdel/(nx-0.99999d0) or $
         min(yy(*,1:ny-1)-yy(*,0:ny-2)) lt ydel/(ny-0.99999d0) then $
           irregular = 1

if not irregular then begin
   ureg=u
   vreg=v
endif else begin
   ; Check if the triangulation has already been done
   newx=1
   if keyword_set(triangles) and $
       n_elements(xxold) eq nx*ny and n_elements(yyold) eq nx*ny then $
          if max(abs(xx-xxold))+max(abs(yy-yyold)) eq 0 then newx=0
   if newx then begin
      print,'Triangulating in VECTOR ...'
      triangulate,float(xx),float(yy),triangles
      xxold=xx
      yyold=yy
   endif
   ; Calculate optimal size for regular grid. Overwrite nx and ny!
   ; (i) Number of regular cells is the same as number of irregular points
   ; (ii) cells should be squares
   aspect = xdel/ydel
   nx     = floor(sqrt(aspect*nx*ny))
   ny     = floor(nx/aspect)

   ; Interpolate velocity arrays onto a REGULAR grid with nx'*ny' spacing
   dx     = xdel/(nx-1.00001d0)
   dy     = ydel/(ny-1.00001d0)

   ; Interpolate
   ureg=trigrid(xx,yy,u,triangles,[dx,dy],[xmin,ymin,xmax,ymax])
   vreg=trigrid(xx,yy,v,triangles,[dx,dy],[xmin,ymin,xmax,ymax])
endelse

minval=maxval*1e-4

; Calculate the coordinates forming the arrows: X(nvecs,nsteps+4,2)
; --> nvecs arrows to be drawn as 
; --> nsteps straight line segments that trace out the vector field
;     with 3 additional line segments forming the arrow head. 
;     The nsteps+3 segments are defined by altogether nsteps+4 points.
; --> third dimension is 0 and 1 for X and Y

; Starting positions
;
x=fltarr(nvecs,nsteps+4,2)
x(*,0,*)=x0

; The time step based on the length of the longest arrow and the maxval speed
;
dt=length*sqrt(xdel^2+ydel^2)/nsteps/maxval

; Integrate the velocity field for nsteps
;
for i=1,nsteps do begin
   xt=(nx-1)*(x(*,i-1,0)-xmin)/xdel
   yt=(ny-1)*(x(*,i-1,1)-ymin)/ydel
   ut=interpolate(ureg,xt,yt)
   vt=interpolate(vreg,xt,yt)
   if irregular and i eq 1 then for ivec=0,nvecs-1 do begin
      if random and abs(ut(ivec)) lt minval and abs(vt(ivec)) lt minval then $
      begin
         xivec=randomu(seed)         & yivec=randomu(seed)
         x(ivec,0,0)=xmin+xdel*xivec & x(ivec,0,1)=ymin+ydel*yivec
         xivec=xivec*(nx-1)          & yivec=yivec*(ny-1)
         ut(ivec)=interpolate(ureg,xivec,yivec)
         vt(ivec)=interpolate(vreg,xivec,yivec)
      endif
   endfor

   ; Second order for stream line integration
   if nsteps gt 5 then begin
      xt1=xt+0.5*dt*ut*(nx-1)/xdel
      yt1=yt+0.5*dt*vt*(ny-1)/ydel
      ut=interpolate(ureg,xt1,yt1)
      vt=interpolate(vreg,xt1,yt1)
   endif

   x(*,i,0)=x(*,i-1,0)+ut*dt
   x(*,i,1)=x(*,i-1,1)+vt*dt
endfor

; Put the position of the arrowhead into the elements nsteps+1..nsteps+3.
; The arrow head consists of two "wings" starting from the tip of the arrow.
; The wings are at a 30 degree angle (taking into account the aspect ratio).
; The projected length of the wings onto the arrow is head*nsteps*last_segment.
; The last segment is between x(*,nsteps,*) and x(*,nsteps-1,*).

aspect = xdel/ydel
tant   = tan(30/!radeg)
scal   = head*nsteps

i = nsteps   &   h = nsteps-1

dx = scal*(x(*,i,0)-x(*,h,0))	   ; Vector of the projected wings
dy = scal*(x(*,i,1)-x(*,h,1))
xm = x(*,i,0)-dx                   ; The coordinates of the midpoint
ym = x(*,i,1)-dy                   ; between the two 'wings'
dx = dx*tant/aspect
dy = dy*tant*aspect                ; Contract by tan(theta) and use aspect.

x(*,i+1,0) = xm-dy                 ; (dx,dy) is rotated to (-dy,dx) and added
x(*,i+2,0) = x(*,i,0)		   ; to and subtracted from (xm,ym).
x(*,i+3,0) = xm+dy

x(*,i+1,1) = ym+dx
x(*,i+2,1) = x(*,i,1)
x(*,i+3,1) = ym-dx

; Put new arrow positions into X0 for next call based on dynamic [0..nsteps]
x0(*,*) = x(*,dynamic,*)

; Draw box
plot,[xmin,xmax],[ymin,ymax],/nodata,XSTYLE=1,YSTYLE=1,NOERASE=noerase

; Draw arrows
if keyword_set(white) then begin
    tvlct,bytarr(256,3)+255       ; change to white
    for i=0,nvecs-1 do plots,x(i,*,0),x(i,*,1),NOCLIP=0,color=127
    common colors                 ; restore colors
    tvlct,r_curr,g_curr,b_curr
endif else $
  for i=0,nvecs-1 do plots,x(i,*,0),x(i,*,1),NOCLIP=0

return
end
;=============================================================================

