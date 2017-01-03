; This code is a copyright protected software (c) 2002- University of Michigan 
;
; Definitions and/or default values for global variables

; System variables that can get corrupted if an animation is interrupted
!x.tickname=strarr(60)
!y.tickname=strarr(60)

; Confirmation for set parameters
common ask_param, $
   doask
doask=0

common fits_param, $
   noresize
noresize=0  ; Keep original size of fits image

; Parameters for .r getpict
common getpict_param, $
   filename, nfile, filenames, filetypes, npictinfiles, npict
filename=''         ; space separated string of filenames. May contain *, []
nfile=0             ; number of files
filenames=0         ; array of filenames
filetype=''         ; file types (binary, real4, ascii...)
npictinfiles=0      ; number of pictures in each file
npict=0             ; index of snapshot to be read

; Parameters for .r plotfunc
common plotfunc_param, $
   func, nfunc, funcs, funcs1, funcs2, plotmode, plotmodes, nplot, $
   plottitle, plottitles, plottitles_file, $
   timetitle, timetitleunit, timetitlestart, $
   autorange, autoranges, noautorange, fmin, fmax, $
   axistype, bottomline, headerline
func=''             ; space separated list of functions to be plotted
nfunc=0             ; number of functions to be plotted
funcs=''            ; array of function names
funcs1=''           ; array of first  components of vectors functions
funcs2=''           ; array of second components of vectors functions
plotmode='plot'     ; space separated list of plot modes
plotmodes=''        ; array of plot modes
nplot=0             ; number of subplots (overplot functions count as 1)
plottitle='default' ; semicolon separated list of titles
plottitles=''       ; array of plot titles
plottitles_file=''  ; array of plottitle strings per file
timetitle=''        ; set to format string to plot time as title for time series
timetitleunit=0     ; set to number of seconds in time unit
timetitlestart=0    ; set to initial time to be subtracted (in above units)
autorange='y'       ; function ranges fmin and fmax set automatically or by hand
autoranges=''       ; array of autorange values
noautorange=0       ; true if all autoranges are 'n'
fmin=0              ; array of minimum values for each function
fmax=0              ; array of maximum values for each function
axistype='coord'    ; 'cells' or 'coord'
headerline=0        ; Number of items to show at the top
bottomline=3        ; Number of items or a string to show at the bottom

; Animation parameters for the movie
common animate_param, $
   firstpict, dpict, npictmax, savemovie, wsubtract, timediff
firstpict=1   ; a scalar or array (per file) of the index of first frame
dpict=1       ; a scalar or array (per file) of distance between frames
npictmax=500  ; maximum number of frames in an animation
savemovie='n' ; save movie into 'ps', 'png', 'tiff', 'bmp', 'jpeg' files
wsubtract=0   ; Array subtracted from w during animation
timediff=0    ; take time derivative of w during animation if timediff=1

; Parameters for .r slice
common slice_param, $
   firstslice, dslice, nslicemax, slicedir, dyslicelabel, $
   x2d, y2d, var2d, grid2d
firstslice=1      ; index of first slice
dslice=1          ; stride between slices
nslicemax=500     ; maximum number of slices shown
slicedir=0        ; 
dyslicelabel=0.98 ; position of bottom label (?)
x2d=0             ; 1st coordinate of the slice
y2d=0             ; 2nd coordinate of the slice 
var2d=0           ; variable names for the slice
grid2d=0          ; grid indexes for the slice

; Transformation parameters for irregular grids 
common transform_param, $
   usereg, dotransform, transform, nxreg, xreglimits, wregpad, $
   nxreg_old, xreglimits_old, triangles, $
   symmtri
usereg=0          ; use wreg and xreg instead of w and x
dotransform='n'   ; do transform with .r plotfunc?
transform='n'     ; transformation 'none', 'regular', 'my', 'polar', 'unpolar'
nxreg=[0,0]       ; size of transformed grid
xreglimits=0      ; limits of transformed grid [xmin, ymin, xmax, ymax]
wregpad=0         ; array of values used in "padding" the regular arrays
nxreg_old=0       ; previous transformation grid
xreglimits_old=0  ; previous limits of transform grid
triangles=0       ; triangulation saved from previous transform
symmtri=0         ; use symmetric triangulation during transformation?

common vector_param, $
   nvector, vectors
nvector=0         ; number of vector variables
vectors=0         ; index of first components of vector variables

; Parameters for getlog
common getlog_param, $
   logfilename, logfilenames
logfilename=''    ; space separated string of filenames. May contain *, []
logfilenames=0    ; array of log filenames

common log_data, $
   timeunit, $
   wlog , logtime , wlognames , $
   wlog1, logtime1, wlognames1, $
   wlog2, logtime2, wlognames2, $
   wlog3, logtime3, wlognames3, $
   wlog4, logtime4, wlognames4, $
   wlog5, logtime5, wlognames5, $
   wlog6, logtime6, wlognames6, $
   wlog7, logtime7, wlognames7, $
   wlog8, logtime8, wlognames8, $
   wlog9, logtime9, wlognames9

wlog=0       ; data array from logfile
logtime=0    ; time array from logfile
wlognames='' ; array of log data names
timeunit='h' ; set to '1' (unitless), 's' (second), 'm' (minute), 'h' (hour) 
	     ;        'millisec', 'microsec', 'ns' (nanosec)

; Parameters for plotlog
common plotlog_param, $
   log_spacex,log_spacey, logfunc, title, xtitle, ytitles, $
   xrange, yranges, colors, linestyles, symbols, smooths, dofft, $
   legends, legendpos

log_spacex=5 ; horizontal distance between log plots measured in character size
log_spacey=5 ; vertical distance between log plots measured in character size
logfunc=''   ; space separated list of log variables in wlogname(s)
title=0      ; set to a string with the title
xtitle=0     ; set to a string with the time title
ytitles=0    ; set to a string array with the function names
xrange=0     ; set to a [min,max] array for the time range
yranges=0    ; set to a [[min1,max1], [min2,max2] ...] for function ranges
colors=255   ; set to an array with colors for each function
linestyles=0 ; set to an array with line styles for each function
symbols=0    ; set to an array with symbols for each function
smooths=0    ; set to an array with smoothing width for each logfile
dofft=0      ; set to 1 to do an FFT transform on the functions
legends=''   ; legends for the lines, defaults are logfilenames
legendpos=0  ; position for the legends: [xmin, xmax, ymin, ymax]

; arrays containing plot data
common plot_data, $
   grid, $
   x, w, xreg, wreg, $
   x0,w0, x1,w1, x2,w2, x3,w3, x4,w4, x5,w5, x6,w6, x7,w7, x8,w8, x9,w9, $
   wreg0, wreg1, wreg2, wreg3, wreg4, wreg5, wreg6, wreg7, wreg8, wreg9
grid=0              ; index array for the whole grid to be used for cuts
x=0                 ; coordinate array of the last read snapshot
w=0                 ; data array of the last read snapshot
xreg=0              ; regular grid coordinates
wreg=0              ; regular grid data

; parameters passed to plot_func through common blocks
common plot_param, $
   multiplot, multix, multiy, multidir, plotix, plotiy, $
   plot_spacex, plot_spacey, $
   fixaspect, noerase, $ 
   cut, cut0, plotdim, rcut, rbody, $
   velvector, velpos, velpos0, velrandom, velspeed, velx, vely, veltri, $
   ax, az, contourlevel, linestyle

; multiplot=0 gives the default number of subplots depending on nfile,nfunc
; multiplot=[3,2,0] defines 3 by 2 subplots filled up in vertical order
; multiplot=[2,4,1] defines 2 by 4 subplots filled up in horizontal order
multiplot=0   ;
multix=0      ; number of subplots horizontally
multiy=0      ; number of subplots vertically
multidir=0    ; subplot ordering: horizontal (0) or vertical (1) first  
plotix=0      ; horizontal index of subplot
plotiy=0      ; vertical index of subplot
plot_spacex=3 ; horizontal distance between plots measured in character size
plot_spacey=3 ; vertical distance between plots measured in character size
fixaspect=1   ; Aspect ratio: 0 - not fixed (fill screen), 1 - use coordinates
noerase=0     ; Do not erase before new plot
cut=0         ; index array for the cut
cut0=0        ; index array for the cut without degenerate indices
plotdim=2     ; plot dimensionality after cut is applied
rcut=-1.0     ; radius of cutting out inner part
rbody= -1.    ; radius of the inner body shown as a black circle
velvector=200 ; number of vectors/stream lines per plot
velpos   =0   ; 2 x velvector array with start positions 
velpos0  =0   ; previous velpos
velrandom=0   ; if 1 use random start positions for each frame of animation
velspeed =5   ; speed of moving vectors during animation
velx=0        ; storage for x coordinate of vector positions
vely=0        ; storage for y coordinate of vector positions
veltri=0      ; storage for triangulation
ax=30         ; view angle relative to x axis for surface and shade_surf
az=30         ; view angle relative to z axis for surface and shade_surf
contourlevel=30 ; Number of contour levels for contour and contfill
linestyle=0   ; line style

; store plot function values from plotting and animations
; calculate running max or mean of functions during animation
common plot_store, $
   f, f1, f2, $
   nplotstore, iplotstore, nfilestore, ifilestore, plotstore, timestore
f=0            ; last plot function magnitude
f1=0           ; first  component of last vector plot function
f2=0           ; second component of last vector plot function
nplotstore = 0 ; number of plots stored
iplotstore = 0 ; current plot index
nfilestore = 1 ; number of files stored
ifilestore = 0 ; current file index
plotstore  = 0 ; array of stored data
timestore  = 0 ; array of stored times

; Some useful constants in SI units
common phys_const, kbSI, mpSI, mu0SI, eSI, ReSI, RsSI, AuSI, cSI, e0SI

kbSI   = 1.3807d-23      ; Boltzmann constant
mpSI   = 1.6726d-27      ; proton mass
mu0SI  = 4*!dpi*1d-7     ; vacuum permeability
eSI    = 1.602d-19       ; elementary charge
ReSI   = 6378d3          ; radius of Earth
RsSI   = 6.96d8          ; radius of Sun
AuSI   = 1.4959787d11    ; astronomical unit
cSI    = 2.9979d8        ; speed of light
e0SI   = 1/(mu0SI*cSI^2) ; vacuum permettivity 

; Physical unit names and values in SI units
common phys_units, $
   fixunits, typeunit, xSI, tSI, rhoSI, uSI, pSI, bSI, jSI, Mi, Me, $
   Qi, Qe, gamma, gammae, clight

fixunits   = 0             ; fix units (do not overwrite) if true
typeunit   = 'NORMALIZED'  ; 'SI', 'NORMALIZED', 'PIC', 'PLANETARY', 'SOLAR'
xSI        = 1.0           ; distance unit in SI
tSI        = 1.0           ; time unit in SI
rhoSI      = 1.0           ; density unit in SI
uSI        = 1.0           ; velocity unit in SI
pSI        = 1.0           ; pressure unit in SI
bSI        = sqrt(mu0SI)   ; magnetic unit in SI
jSI        = 1/sqrt(mu0SI) ; current unit in SI
Mi         = 1.0           ; Ion mass in amu
Me         = 1/1836.15     ; Electron mass in amu
Qi         = 1.0           ; Ion charge in unit charge
Qe         = -1.0          ; Electron change in unit charge
gamma      = 5./3.         ; Adiabatic index for first fluid
gammae     = 5./3.         ; Adiabatic index for electrons
clight     = cSI           ; (Reduced) speed of light

; conversion factors that are useful to calculate various derived quantities
common phys_convert, $
   ti0, cs0, mu0A, mu0, c0, uH0, op0, oc0, rg0, di0, ld0

ti0  = 1.0 ; ion temperature = ti0*p/rho*Mion
cs0  = 1.0 ; sound speed     = sqrt(cs0*gamma*p/rho)
c0   = 1.0 ; speed of light  = c0
mu0A = 1.0 ; Alfven speed    = sqrt(bb/mu0A/rho)
mu0  = 1.0 ; plasma beta     = p/(bb/2/mu0)
uH0  = 1.0 ; Hall velocity   = uH0*j/rho*Mion
op0  = 1.0 ; plasma freq.    = op0*sqrt(rho)/Mion
oc0  = 1.0 ; cyclotron freq. = oc0*b/mIon
rg0  = 1.0 ; ion gyro radius = rg0*sqrt(p/rho)/b*sqrt(Mion)
di0  = 1.0 ; inertial length = di0*sqrt(rho)*Mion
ld0  = 1.0 ; Debye length    = ld0*sqrt(p)/rho*Mion

; information obtained from the last file header
common file_head, $
   ndim, headline, it, time, gencoord, neqpar, nw, nx, eqpar, $
   variables, wnames

ndim     = 1  ; number of spatial dimensions
headline = '' ; first line often containing physcial unit names
it       = 0  ; time step
time     = 0. ; simulation time
gencoord = 0  ; true for unstructured/non-Cartesian grids
neqpar   = 0  ; number of scalar parameters
nw       = 1  ; number of variables
nx       = 0  ; number of grid cells
eqpar    = 0. ; values of scalar parameters
variables= '' ; array of coordinate/variable/param names
wnames   = '' ; array of variables names

end
