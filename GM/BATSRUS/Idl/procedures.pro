;  Copyright (C) 2002 Regents of the University of Michigan, 
;  portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;
; Written by G. Toth for the Versatile Advection Code and BATSRUS/SWMF
; Some improvements by Aaron Ridley.
;
; Main procedures to (re)set defaults and to read, plot, animate and slice data
;   set_default_values
;   read_data, plot_data, animate_data, slice_data
;   read_log_data, plot_log_data
;
; Procedures for
;
; fixing things if animate or slice crashes
;   reset_axis, slice_data_restore
; reading ascii and binary data produced by VAC, BATSRUS, PWOM, IPIC3D etc:
;    open_file, get_file_types, get_file_head, get_pict, 
;    get_pict_asc, get_pict_bin, get_pict_log, get_log
; showing / overwriting information read from last file:
;    show_head, show_units, set_units
; saving ascii and binary data in the same format as used for input:
;    save_pict, save_log
; reading numbers and strings from input:
;    asknum, askstr, string_to_array, arr2arr, read_plot_param, read_limits
; transforming initial data:
;    read_transform_param, do_transform, do_my_transform, 
;    make_regular_grid, make_polar_grid, make_unpolar_grid,
;    make_sphere_grid, 
;    getaxes
;    interpol_logfiles, interpol_log
; calculating functions of the data
;    get_func, get_limits
; plotting
;    plot_func, plot_grid, plot_log
; calculating cell corners and cell volumes for general 2D grids
;    gengrid
; comparing two w,x or wlog arrays for relative differences
;    compare, rms_logfiles, rel_error
; checking index ranges for functions quadruplet and triplet
;    checkdim
; procedure "quit" as an alias for "exit"
;    quit
;
; Functions for
;
; calculating first derivative in 1D
;    diff1
; calculating derivatives in 2D for Cartesian grids to 2nd,3rd,4th order
;    diff2,diff3,diff4
; calculating the laplace on 1D or 2D Cartesian grids
;    laplace
; calculating minmod limited slope
;    minmod
; calculating symmetric differences with respect to some mirror plane
;    symmdiff
; calculating derivatives in 2D for general grids
;    grad_2d,div,curl,grad_rz,div_rz,curl_rz, filledge,intedge,intedge_rz
; taking a part of an array or coarsen an array
;    triplet, quadruplet, coarsen
; eliminating degenerate dimensions from an array
;    reform2
; converting logfile time or date+time into hours, getting other functions
;    log_time, log_func


;===========================================================================
pro set_default_values

; Definitions and default values for variables in common blocks

; System variables that can get corrupted if an animation is interrupted
  !x.tickname=strarr(60)
  !y.tickname=strarr(60)

; Confirmation for set parameters
  common ask_param, $
     doask
  doask=0
  
; behavior on error: 0: stop in the unit (for debug), 2=return to main
  common debug_param, $
     onerror
  onerror=2

  common fits_param, $
     noresize
  noresize=0                    ; Keep original size of fits image

; Parameters for read_data
  common getpict_param, $
     filename, nfile, filenames, filetypes, npictinfiles, npict
  filename=''          ; space separated list of filenames. May contain *, []
  nfile=0              ; number of files
  filenames=0          ; array of filenames
  filetype=''          ; file types (binary, real4, ascii...)
  npictinfiles=0       ; number of pictures in each file
  npict=0              ; index of snapshot to be read

; Parameters for plot_data
  common plotfunc_param, $
     func, nfunc, funcs, funcs1, funcs2, plotmode, plotmodes, nplot, $
     plottitle, plottitles, plottitles_file, $
     timetitle, timetitleunit, timetitlestart, $
     autorange, autoranges, noautorange, fmin, fmax, $
     axistype, bottomline, headerline
  func=''              ; space separated list of functions to be plotted
  nfunc=0              ; number of functions to be plotted
  funcs=''             ; array of function names
  funcs1=''            ; array of first  components of vectors functions
  funcs2=''            ; array of second components of vectors functions
  plotmode='plot'      ; space separated list of plot modes
  plotmodes=''         ; array of plot modes
  nplot=0              ; number of subplots (overplot functions count as 1)
  plottitle='default'  ; semicolon separated list of titles
  plottitles=''        ; array of plot titles
  plottitles_file=''   ; array of plottitle strings per file
  timetitle=''     ; set to format string to plot time as title for time series
  timetitleunit=0      ; set to number of seconds in time unit
  timetitlestart=0     ; set to initial time to be subtracted (in above units)
  autorange='y'        ; function ranges set automatically or by fmin/fmax
  autoranges=''        ; array of autorange values
  noautorange=0        ; true if all autoranges are 'n'
  fmin=0               ; array of minimum values for each function
  fmax=0               ; array of maximum values for each function
  axistype='coord'     ; 'cells' or 'coord'
  headerline=0         ; Number of items to show at the top
  bottomline=3         ; Number of items or a string to show at the bottom

; Animation parameters for the movie
  common animate_param, $
     firstpict, dpict, npictmax, savemovie, wsubtract, timediff, $
     videosave, videofile, videorate, videoobject, videostream, videotime
  firstpict=1        ; a scalar or array (per file) of the index of first frame
  dpict=1            ; a scalar or array (per file) of distance between frames
  npictmax=500       ; maximum number of frames in an animation
  savemovie='n'      ; save animation frames into ps/png/tiff/bmp/jpeg files
                     ; or into a 'mov/mp4/avi' video file.
  wsubtract=0        ; Array subtracted from w during animation
  timediff=0         ; take time derivative of w during animation if timediff=1
  videosave=0        ; save video?
  videofile='movie'  ; name of video file with extension .savemovie
  videorate=10       ; number of frames per second
  videoobject=0      ; video object
  videostream=0      ; video stream object
  videotime=0        ; length of video

; Parameters for .r slice
  common slice_param, $
     firstslice, dslice, nslicemax, slicedir, dyslicelabel, $
     x3d, w3d, var3d, rbody3d, grid2d
  firstslice=1                  ; index of first slice
  dslice=1                      ; stride between slices
  nslicemax=500                 ; maximum number of slices shown
  slicedir=0                    ; 
  dyslicelabel=0.98             ; position of bottom label (?)
  x3d=0                         ; 3D coordinates saved
  w3d=0                         ; 3D state saved
  var3d=0                       ; 3D variable names saved
  rbody3d=0.0                   ; 3D rBody value saved
  grid2d=0                      ; grid indexes for the slice

; Transformation parameters for irregular grids 
  common transform_param, $
     usereg, dotransform, transform, nxreg, xreglimits, wregpad, $
     nxreg_old, xreglimits_old, triangles, $
     symmtri
  usereg=0         ; use wreg and xreg instead of w and x
  dotransform='n'  ; do transform with plot_data?
  transform='n'    ; transformation 'none', 'regular', 'my', 'polar', 'unpolar'
  nxreg=[0,0]      ; size of transformed grid
  xreglimits=0     ; transformed grid limits [xmin, ymin, xmax, ymax]
  wregpad=0        ; array of values used in "padding" the regular arrays
  nxreg_old=0      ; previous transformation grid
  xreglimits_old=0 ; previous limits of transform grid
  triangles=0      ; triangulation saved from previous transform
  symmtri=0        ; use symmetric triangulation during transformation?

  common vector_param, $
     nvector, vectors
  nvector=0                     ; number of vector variables
  vectors=0                     ; index of first components of vector variables

; Parameters for read_log_data
  common getlog_param, $
     logfilename, logfilenames
  logfilename=''       ; space separated string of filenames. May contain *, []
  logfilenames=0                ; array of log filenames

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

  wlog=0                        ; data array from logfile
  logtime=0                     ; time array from logfile
  wlognames=''                  ; array of log data names
  timeunit='h' ; set to '1' (unitless), 's' (second), 'm' (minute), 'h' (hour) 
                                ;        'millisec', 'microsec', 'ns' (nanosec)

; Parameters for plot_log_data
  common plotlog_param, $
     log_spacex,log_spacey, logfunc, title, xtitle, ytitles, $
     xrange, yranges, timeshifts, $
     colors, linestyles, symbols, smooths, dofft, legends, legendpos

  log_spacex=5 ; horizontal distance around log plots (in character size)
  log_spacey=5 ; vertical distance between log plots (in character size)
  logfunc=''   ; space separated list of log variables in wlogname(s)
  title=''     ; set to a string with the title
  xtitle=0     ; set to a string with the time title
  ytitles=0    ; set to a string array with the function names
  xrange=0     ; set to a [min,max] array for the time range
  yranges=0    ; set to a [[min1,max1], [min2,max2] ...] for function ranges
  timeshifts=0 ; set to array of time shifts per logfile
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
  grid=0                        ; index array to be used for cuts
  x=0                           ; coordinate array of the last read snapshot
  w=0                           ; data array of the last read snapshot
  xreg=0                        ; regular grid coordinates
  wreg=0                        ; regular grid data

; parameters passed to plot_func through common blocks
  common plot_param, $
     multiplot, multix, multiy, multidir, plotix, plotiy, $
     plot_spacex, plot_spacey, showxaxis, showyaxis, showxtitle, showytitle, $
     fixaspect, noerase, $ 
     cut, cut0, plotdim, rcut, rbody, $
     velvector, velpos, velpos0, velrandom, velspeed, velx, vely, veltri, $
     viewanglex, viewanglez, colorlevel, contourlevel, linestyle

; multiplot=0         - default subplot arrangement based on nfile,nfunc
; multiplot=-1        - arrange subplots vertically
; multiplot=3         - 3 subplots horizontally
; multiplot=[3,2,1]   - 3x2 subplots filled up in vertical order
; multiplot=[2,4,0]   - 2x4 subplots filled up in horizontal order
; multiplot=[2,4,0,4] - start with subplot 4 (3rd row left)

  multiplot=0                   ;
  multix=0                      ; number of subplots horizontally
  multiy=0                      ; number of subplots vertically
  multidir=0                    ; subplot order: horizontal (0) or vertical (1)
  plotix=0                      ; horizontal index of subplot
  plotiy=0                      ; vertical index of subplot
  plot_spacex=3                 ; horizontal subplot distance (character size)
  plot_spacey=3                 ; vertical subplot distance (character size)
  showxtitle=0                  ; show x title and axis in all subplots
  showytitle=0                  ; show y title and axis in all subplots
  showxaxis=0                   ; show x axis in all subplots
  showyaxis=0                   ; show y axis in all subplots
  fixaspect=1                   ; fix aspect ratio according to coordinates
  noerase=0                     ; Do not erase before new plot
  cut=0                         ; index array for the cut
  cut0=0                        ; cut array without degenerate indices
  plotdim=2                     ; plot dimensionality after cut is applied
  rcut=-1.0                     ; radius of cutting out inner part
  rbody= -1.                    ; radius of inner body shown as a black circle
  velvector=200                 ; number of vectors/stream lines per plot
  velpos   =0                   ; 2 x velvector array with start positions 
  velpos0  =0                   ; previous velpos
  velrandom=0                   ; use new random start positions in animation
  velspeed =5                   ; speed of moving vectors during animation
  velx=0                        ; storage for x coordinate of vector positions
  vely=0                        ; storage for y coordinate of vector positions
  veltri=0                      ; storage for triangulation
  viewanglex=30                 ; view angle relative to x axis for surface/shade_surf
  viewanglez=30                 ; view angle relative to z axis for surface/shade_surf
  colorlevel=30                 ; Number of color levels for contfill/contbar
  contourlevel=30               ; Number of contour levels for contour
  linestyle=0                   ; line style for plot

; store plot function values from plotting and animations
; calculate running max or mean of functions during animation
  common plot_store, $
     f, f1, f2, $
     nplotstore, iplotstore, nfilestore, ifilestore, plotstore, timestore
  f=0                           ; last plot function magnitude
  f1=0                          ; first  component of last vector plot function
  f2=0                          ; second component of last vector plot function
  nplotstore = 0                ; number of plots stored
  iplotstore = 0                ; current plot index
  nfilestore = 1                ; number of files stored
  ifilestore = 0                ; current file index
  plotstore  = 0                ; array of stored data
  timestore  = 0                ; array of stored times

; Some useful constants in SI units
  common phys_const, kbSI, mpSI, mu0SI, eSI, ReSI, RsSI, AuSI, cSI, e0SI

  kbSI   = 1.3807d-23           ; Boltzmann constant
  mpSI   = 1.6726d-27           ; proton mass
  mu0SI  = 4*!dpi*1d-7          ; vacuum permeability
  eSI    = 1.602d-19            ; elementary charge
  ReSI   = 6378d3               ; radius of Earth
  RsSI   = 6.96d8               ; radius of Sun
  AuSI   = 1.4959787d11         ; astronomical unit
  cSI    = 2.9979d8             ; speed of light
  e0SI   = 1/(mu0SI*cSI^2)      ; vacuum permettivity 

; Physical unit names and values in SI units
  common phys_units, $
     fixunits, typeunit, xSI, tSI, rhoSI, uSI, pSI, bSI, jSI, Mi, Me, $
     Qi, Qe, gamma, gammae, clight

  fixunits   = 0                ; fix units (do not overwrite) if true
  typeunit   = 'NORMALIZED'     ; 'SI'/'NORMALIZED'/'PIC'/'PLANETARY'/'SOLAR'
  xSI        = 1.0              ; distance unit in SI
  tSI        = 1.0              ; time unit in SI
  rhoSI      = 1.0              ; density unit in SI
  uSI        = 1.0              ; velocity unit in SI
  pSI        = 1.0              ; pressure unit in SI
  bSI        = sqrt(mu0SI)      ; magnetic unit in SI
  jSI        = 1/sqrt(mu0SI)    ; current unit in SI
  Mi         = 1.0              ; Ion mass in amu
  Me         = 1/1836.15        ; Electron mass in amu
  Qi         = 1.0              ; Ion charge in unit charge
  Qe         = -1.0             ; Electron change in unit charge
  gamma      = 5./3.            ; Adiabatic index for first fluid
  gammae     = 5./3.            ; Adiabatic index for electrons
  clight     = cSI              ; (Reduced) speed of light

; conversion factors that are useful to calculate various derived quantities
  common phys_convert, $
     ti0, cs0, mu0A, mu0, c0, uH0, op0, oc0, rg0, di0, ld0, vec0

  ti0  = 1.0                    ; ion temperature = ti0*p/rho*Mi
  cs0  = 1.0                    ; sound speed     = sqrt(cs0*gamma*p/rho)
  c0   = 1.0                    ; speed of light  = c0
  mu0A = 1.0                    ; Alfven speed    = sqrt(bb/mu0A/rho)
  mu0  = 1.0                    ; plasma beta     = p/(bb/2/mu0)
  uH0  = 1.0                    ; Hall velocity   = uH0*j/rho*Mi
  op0  = 1.0                    ; plasma freq.    = op0*sqrt(rho)/Mi
  oc0  = 1.0                    ; cyclotron freq. = oc0*b/mIon
  rg0  = 1.0                    ; ion gyro radius = rg0*sqrt(p/rho)/b*sqrt(Mi)
  di0  = 1.0                    ; inertial length = di0*sqrt(rho)*Mi
  ld0  = 1.0                    ; Debye length    = ld0*sqrt(p)/rho*Mi
  vec0 = [0.0, 0.0, 1.0]        ; vector to define field aligned coordinates

; information obtained from the last file header
  common file_head, $
     ndim, headline, it, time, gencoord, neqpar, nw, nx, eqpar, $
     variables, wnames

  ndim     = 1                  ; number of spatial dimensions
  headline = ''                 ; 1st line often containing physcial unit names
  it       = 0                  ; time step
  time     = 0.                 ; simulation time
  gencoord = 0                  ; true for unstructured/non-Cartesian grids
  neqpar   = 0                  ; number of scalar parameters
  nw       = 1                  ; number of variables
  nx       = 0                  ; number of grid cells
  eqpar    = 0.                 ; values of scalar parameters
  variables= ''                 ; array of coordinate/variable/param names
  wnames   = ''                 ; array of variables names

; The byte arrays in the colors common block are set by loadct
  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

end
;===========================================================================
pro read_data
;
;    Read the npict-th snapshot from an ascii or binary data file into
;    the x (coordinates) and w (data) arrays. 
;    If dotransfrom='y' the data is transformed according to 
;    the transformation parameters, potentially into xreg, wreg.
;
;    Usage: 
;
; filename='...' ; set file to read from (optional)
; npict=...      ; set snapshot index (optional)
; read_data
;
;    read_data will prompt you for "filename(s)" and "npict"
;    unless they are already set. Previous settings can be erased by 
;
; set_defaults
;
;    or modified explicitly, e.g.:
;
; filename='data/example.ini'
; npict=1
;
;    The "x" and "w" arrays and the header info will be read from the file. 
;
;    If a file is read with generalized coordinates, "gencoord=1" is set,
;    and the original data is transformed according to the "transform"
;    string variable into "xreg" and "wreg".
;
;    The same npict-th snapshot can be read from up to 10 files by e.g. setting
;
; filename='data/file1.ini data/file2.out'
;
;    In this case the data is read into x0,w0 and x1,w1 for the two files,
;    and possibly transformeed into wreg0,wreg1.
;
;    To plot a variable with IDL functions, type e.g.:
;
; plot, x(*,0), w(*,5,2)
;
;    or use the procedure
;
; plot_data
;
;===========================================================================

  common debug_param & on_error, onerror

  common ask_param, doask
  common getpict_param
  common plot_data
  common file_head
  common transform_param

  nfile=0
  askstr,'filename(s)   ',filename, doask
  string_to_array, filename, filenames, nfile, /wildcard

  if not keyword_set(filenames) then begin
     print,'Error in read_data: no matching filename was found.'
     retall
  endif

  if nfile gt 10 then begin
     print,'Error in read_data: cannot handle more than 10 files.'
     print,'nfile     = ', nfile
     print,'filenames = ', filenames
     retall
  endif
  get_file_types
  print,'filetype(s)   =','',filetypes
  print,'npictinfile(s)=',npictinfiles
  if max(npictinfiles) eq 1 then npict=1
  asknum,'npict',npict,doask
  print

  for ifile=0,nfile-1 do begin

     ;; Read data from file ifile

     print

     open_file,10,filenames(ifile),filetypes(ifile)
     get_pict,10,filenames(ifile),filetypes(ifile),npict<npictinfiles(ifile),$
              error

     show_head, ifile

     if nfile gt 1 then begin
        case ifile of
           0: begin
              w0=w
              x0=x
           end
           1: begin
              w1=w
              x1=x
           end
           2: begin
              w2=w
              x2=x
           end
           3: begin
              w3=w
              x3=x
           end
           4: begin
              w4=w
              x4=x
           end
           5: begin
              w5=w
              x5=x
           end
           6: begin
              w6=w
              x6=x
           end
           7: begin
              w7=w
              x7=x
           end
           8: begin
              w8=w
              x8=x
           end
           9: begin
              w9=w
              x9=x
           end
        endcase
        print,'Read x',ifile,' and w',ifile,FORMAT='(a,i1,a,i1)'
     endif else print,'Read x and w'

     read_transform_param

     do_transform

     if usereg then begin
        if nfile eq 1 then $
           print,'...transform to xreg and wreg' $
        else print,'...transform to xreg and wreg',ifile,FORMAT='(a,i1)'

        if nfile gt 1 then case ifile of
           0: wreg0 = wreg
           1: wreg1 = wreg
           2: wreg2 = wreg
           3: wreg3 = wreg
           4: wreg4 = wreg
           5: wreg5 = wreg
           6: wreg6 = wreg
           7: wreg7 = wreg
           8: wreg8 = wreg
           9: wreg9 = wreg
        endcase
     endif
  endfor
  close,10

; Produce a wnames from the last file
  wnames=variables(ndim:ndim+nw-1)

end
;===========================================================================
pro plot_data

;    Use the x and w arrays (usually read by read_data or animate_data) and
;    plot one or more functions of w using different plotting routines. 
;    The functions are defined in the "Idl/funcdef.pro" file.
;
;    For generalized coordinates the variables are interpolated from the 
;    irregular grid onto a regular one.
;
;    A subset can be cut from the grid by using the "cut" index array, e.g.:
;    cut=grid(10:30,*), where "grid" contains the full index array.
;    for the regular grid. The grid array is only defined after animate ran.
;
;    Usage:
;
; w=w2 & x=x2    ; set x and w (optional)
; func='....'    ; function(s) to plot (optional)
; plotmode='...' ; plot mode(s) (optional)
; plot_data
;
;    Output can be directed to PS, EPS and/or PDF files:
;
; set_device,'filename.eps'
; plot_data
; close_device,/pdf
;
;===========================================================================

  common debug_param & on_error, onerror

  common getpict_param
  common plot_param
  common plotfunc_param
  common file_head
  common plot_store

  if not keyword_set(nfile) then read_data

  if nfile gt 1 then begin
     print,'More than one files were read...'
     print,'Probably w is from file ',filenames(nfile-1)
     nfile=1
  endif

  print,'======= CURRENT PLOTTING PARAMETERS ================'
  print,'colorlevel=',colorlevel,', contourlevel=',contourlevel,$
        ', velvector=',velvector,', velspeed (0..5)=',velspeed,$
        FORMAT='(a,i3,a,i3,a,i4,a,i2)'
  print,'multiplot=',multiplot
  print,'axistype (coord/cells)=',axistype,', fixaspect= ',fixaspect,$
        FORMAT='(a,a,a,i1)'
  print,'bottomline=',bottomline,', headerline=',headerline,$
        FORMAT='(a,i1,a,i1)'

  if keyword_set(cut) then help,cut
  if keyword_set(velpos) then help,velpos
  velpos0=velpos

                                ; Read plotting and transforming parameters

  print,'======= PLOTTING PARAMETERS ========================='
  print,'wnames                     =',wnames

  read_plot_param

  help,nx

  read_transform_param

  do_transform

  print,'======= DETERMINE PLOTTING RANGES ==================='

  read_limits

  if noautorange eq 0 then begin
     get_limits,1

     print
     for ifunc=0,nfunc-1 do $
        print,'Min and max value for ',funcs(ifunc),':',fmin(ifunc),fmax(ifunc)
  endif

  ;;===== DO PLOTTING IN MULTIX * MULTIY MULTIPLE WINDOWS

  if keyword_set(multiplot) then begin
     if n_elements(multiplot) eq 1 then begin
        if multiplot gt 0 then       !p.multi=[0,multiplot,1 ,0,1] $
        else if multiplot eq -1 then !p.multi=[0,1,nplot     ,0,1] $
        else                         !p.multi=[0,1,-multiplot,0,1]
     endif else if n_elements(multiplot) eq 4 then $
        !p.multi = [multiplot(0)*multiplot(1)-multiplot(3),multiplot(0:1),0,multiplot(2)] $
     else $
        !p.multi=[0,multiplot(0:1),0,multiplot(2)]
     multix=!p.multi(1)
     multiy=!p.multi(2)
  endif else begin
     multix=long(sqrt(nplot-1)+1)
     multiy=long((nplot-1)/multix+1)
     !p.multi=[0,multix,multiy,0,0]
  endelse

  if not noerase then erase

  if velrandom then velpos=0

  if !d.name eq 'X' and !d.window ge 0 then wshow

  plot_func
  
  putbottom,1,1,0,0,bottomline,nx,it,time
  putheader,1,1,0,0,headerline,headline,nx

  print
  !p.multi=0
  !p.title=''

  ;; Restore velpos array
  velpos=velpos0 & velpos0=0

end
;===========================================================================
pro animate_data

;    Written by G. Toth for the Versatile Advection Code.
;
;    Read pictures from one or more ascii or binary data files and 
;    plot or animate one or more functions of w using different plotting 
;    routines. The functions are defined in the "Idl/funcdef.pro" file.
;
;    For generalized coordinates the variables can be interpolated from the 
;    irregular grid onto a regular one, or to polar coordinates.
;
;    A subset can be cut from the grid by using the "cut" index array, e.g.:
;    cut=grid(10:30,*), where "grid" contains the full index array.
;    for the regular grid. The grid array is only defined after animate ran.
;
;    Usage:
; 
; filename='... ; files to read (optional)
; func='...     ; functions to plot (optional)
; animate_data
;
;===========================================================================

  common debug_param & on_error, onerror

  common animate_param
  common getpict_param
  common file_head
  common ask_param
  common plot_param
  common plot_data
  common plotfunc_param
  common colors

  ;; Initialize storage for running maxima and averages
  iplotstore = 0
  plotstore  = 0
  timestore  = 0

  videosave = savemovie eq 'mp4' or savemovie eq 'avi' or savemovie eq 'mov'

  print,'======= CURRENT ANIMATION PARAMETERS ================'
  print,'firstpict=',firstpict,', dpict=',dpict,', npictmax=',npictmax, $
        FORMAT='(a,'+string(n_elements(firstpict))+'i4,a,' $
        +string(n_elements(dpict))+'i4,a,i4)'
  print,'savemovie (n/avi/mp4/mov/ps/png/tiff/bmp/jpeg)=',savemovie
  if videosave then print,'videofile=',videofile,'(.',savemovie,')',$
                          ', videorate=',videorate
  print,'colorlevel=',colorlevel,', contourlevel=',contourlevel,$
        ', velvector=',velvector,', velspeed (0..5)=',velspeed,$
        FORMAT='(a,i3,a,i3,a,i4,a,i2)'
  if keyword_set(multiplot) then begin
     siz=size(multiplot)
     ;; scalar multiplot value is converted to a row (+) or a column (-)
     if siz(0) eq 0 then begin
        if multiplot gt 0 then multiplot=[multiplot,1,1] $
        else                   multiplot=[1,-multiplot,1]
     endif
     print,'multiplot= ',multiplot,', axistype (coord/cells)=',axistype,$
           ', fixaspect= ',fixaspect,$
           FORMAT='(a,"[",i2,",",i2,",",i2,"]",a,a,a,i1)'
  endif else $
     print,'multiplot= 0 (default), axistype (coord/cells)=',axistype,$
           ', fixaspect= ',fixaspect,$
           FORMAT='(a,a,a,i1)'
  print,'bottomline=',bottomline,', headerline=',headerline,$
        FORMAT='(a,i1,a,i1)'

  if keyword_set(cut) then help,cut
  if keyword_set(wsubtract) then help,wsubtract
  if keyword_set(velpos) then help,velpos
  velpos0 = velpos

  print,'======= FILE DESCRIPTION ============================'
  nfile=0
  if filename eq '' or doask then $
     askstr,'filename(s)   ',filename,doask

  string_to_array,filename,filenames,nfile,/wildcard

  get_file_types

  print, 'filenames     =', filenames
  print, 'filetype(s)   =', filetypes
  print, 'npictinfile(s)=', npictinfiles

  ;; Extend firstpict and dpict into arrays of size nfile
  arr2arr,firstpict,nfile
  arr2arr,dpict,nfile

  ;;====== OPEN FILE(S) AND READ AND PRINT HEADER(S)

  anygencoord=0
  for ifile=0,nfile-1 do begin
     open_file,10,filenames(ifile),filetypes(ifile)
     get_file_head,10,filenames(ifile),filetypes(ifile)
     anygencoord=anygencoord or gencoord
     print,         'headline                  =',strtrim(headline,2)
     print,FORMAT='("variables                 =",100(a," "),$)',variables
     print,FORMAT='(" (ndim=",i2,", nw=",i2,")")',ndim,nw
  endfor

  print,'======= PLOTTING PARAMETERS ========================='
  read_plot_param

  read_transform_param

  print,'======= DETERMINE PLOTTING RANGES ==================='

  read_limits

  if noautorange then begin
     npict = min( (npictinfiles-firstpict)/dpict + 1 )
     if npict gt npictmax then npict=npictmax
     if npict lt 0 then npict=0
  endif else begin
     npict=0
     for ifile=0,nfile-1 do $
        open_file,ifile+10,filenames(ifile),filetypes(ifile)
     error=0
     while npict lt npictmax and not error do begin

        for ifile = 0, nfile-1 do begin

           if npict eq 0 then nextpict=firstpict(ifile) $
           else               nextpict=dpict(ifile)

           ;; IPIC3D reader always counts from the beginning
           if filetypes(ifile) eq 'IPIC3D' then $
              nextpict = firstpict(ifile) + npict*dpict(ifile)

           get_pict, ifile+10, filenames(ifile), filetypes(ifile), nextpict, err

           if keyword_set(wsubtract) then w=w-wsubtract

           if keyword_set(timediff) then begin
              if npict eq 0 then begin
                 timeprev = time
                 wprev = w
                 w = 0.0*w
              endif else begin
                 w = (w - wprev)/(time - timeprev)
                 wprev = wprev + w*(time - timeprev)
                 timeprev = time
              endelse
           endif

           wnames=variables(ndim:ndim+nw-1)
           error=err or error

           if not error then begin
              do_transform,ifile

              first= npict eq 0 and ifile eq 0
              get_limits, first

              if ifile eq nfile-1 then begin
                 if npict eq 0 then print,FORMAT='("ipict:    ",$)'
                 npict=npict+1
                 print,FORMAT='(i4,$)',npict
              endif
           endif
        endfor
     endwhile
     print
     for ifunc=0,nfunc-1 do $
        print,'Min and max value for ',funcs(ifunc),':',fmin(ifunc),fmax(ifunc)

  endelse
  print,'npict=',npict
  if npict eq 0 then begin
     print,'There are no frames to animate! Check the following settings:'
     print,'   npictinfiles=',npictinfiles
     print,'   firstpict   =',firstpict
     print,'   dpict       =',dpict
     print,'   npictmax    =',npictmax
     if min(npictinfiles - firstpict) lt 0 then $
        print,'   firstpict is larger than npictinfiles for some files!' 
     retall
  endif

  ;;===== DO ANIMATION IN MULTIX * MULTIY MULTIPLE WINDOWS

  if keyword_set(multiplot) then begin
     multix=multiplot(0)
     multiy=multiplot(1)
     multidir=multiplot(2)
     npict1=(multix*multiy)/(nplot*nfile)
     if npict1 eq 0 then npict1=1
  endif else if nfile eq 1 then begin
     multix=long(sqrt(nplot-1)+1)
     multiy=long((nplot-1)/multix+1)
     multidir=0
     npict1=1
  endif else begin
     multix=nfile
     multiy=nplot
     multidir=1
     npict1=1
  endelse

  if videosave then begin
     videoobject = IDLffVideoWrite(videofile+'.'+savemovie)
     videostream = videoobject.AddVideoStream(!d.x_size,!d.y_size,videorate)
  endif else begin
     if savemovie ne 'n' then spawn,'/bin/mkdir -p Movie'
     if savemovie eq 'ps' then set_plot,'PS',/INTERPOLATE
  endelse

  doanimate= npict gt npict1 and !d.name eq 'X'
  if !d.name eq 'X' then begin
     if !d.window lt 0 then window
     wshow
  endif
  if doanimate then xinteranimate,set=[!d.x_size,!d.y_size,(npict-1)/npict1+1]

  ipict=0
  ipict1=0
  iplot=0
  for ifile=0,nfile-1 do open_file,ifile+10,filenames(ifile),filetypes(ifile)
  error=0
  while ipict lt npict and not error do begin
     if ipict1 eq 0 then begin
        if not keyword_set(noerase) then erase
        !p.multi=[0,multix,multiy,0,multidir]
        if savemovie eq 'ps' then $
           device,filename='Movie/'+string(FORMAT='(i4.4)',iplot+1)+'.ps',$
                  XSIZE=24,YSIZE=18,/LANDSCAPE,/COLOR,BITS=8
     endif

     if ipict eq 0 then print,FORMAT='("ipict:    ",$)'
     print,FORMAT='(i4,$)',ipict+1

     if n_elements(velpos) gt 1 and velrandom gt 0 then begin
        ;; reset a subset of the vector positions to a random value
        ii = velrandom*indgen(velvector/velrandom) + (ipict mod velrandom)
        velpos(ii,0) = randomu(seed,n_elements(ii),/double)*1e6
        velpos(ii,1) = randomu(seed,n_elements(ii),/double)*1e6
     endif

     for ifile=0,nfile-1 do begin

        if npict gt 1 or nfile gt 1 or noautorange then begin

           if ipict eq 0 then nextpict=firstpict(ifile) $
           else               nextpict=dpict(ifile)

           ;; IPIC3D reader always counts from the beginning
           if filetypes(ifile) eq 'IPIC3D' then $
              nextpict = firstpict(ifile) + ipict*dpict(ifile)

           get_pict, ifile+10, filenames(ifile),filetypes(ifile), nextpict, err

           error=error or err
        endif

        if not error then begin

           if keyword_set(wsubtract) then w=w-wsubtract

           if keyword_set(timediff) then begin
              if ipict eq 0 then begin
                 timeprev = time
                 wprev = w
                 w = 0.0*w
              endif else begin
                 w = (w - wprev)/(time - timeprev)
                 wprev = wprev + w*(time - timeprev)
                 timeprev = time
              endelse
           endif

           wnames=variables(ndim:ndim+nw-1)

           do_transform

           linestyle=0
           if multix*multiy lt nplot*nfile then $
              linestyle=fix(nplot*ifile)/(multix*multiy)

           if(keyword_set(timetitle))then begin
              t = time
              if(keyword_set(timetitleunit))then  t = t / timetitleunit
              if(keyword_set(timetitlestart))then t = t - timetitlestart
              plottitles(*) = string(format=timetitle,t)
           endif

           if(keyword_set(plottitles_file))then begin
              plottitle = plottitles_file(ifile)
              string_to_array,plottitle,plottitles,nfunc,';'
           end

           nfilestore = nfile
           ifilestore = ifile

           plot_func

           if npict1 le 1 then begin
              putbottom,multix,multiy,ifile,0,bottomline,nx,it,time
              putheader,nfile,1,ifile,0,headerline,headline,nx
           endif else begin
              putbottom,multix,multiy,plotix,plotiy,bottomline,nx,it,time
              if ipict1 eq 0 then $
                 putheader,nfile,1,ifile,0,headerline,headline,nx
           endelse
        endif

     endfor

     if ipict1 eq npict1-1 or ipict eq npict-1 then begin
        if doanimate then $
           xinteranimate,frame=iplot,window=!d.window
        if videosave then begin
           videotime = videoobject.put(videostream, tvrd(true=1))
        endif else if savemovie eq 'ps' then begin
           print,FORMAT='(" (Movie/",i4.4,".ps)",$)',iplot+1
           device,/close
        endif else if savemovie ne 'n' and !d.name eq 'X' then begin
           imagefile=string(FORMAT='("Movie/",i4.4,".",a)',iplot+1,savemovie)
           print,FORMAT='("(",a,")",$)',imagefile
           write_image, imagefile, savemovie, $
                       tvrd( order=(savemovie eq 'tiff'), true=1), $
                       r_curr, g_curr, b_curr
        endif
     endif

     ipict1=ipict1+1
     if ipict1 ge npict1 then begin
        ipict1=0
        iplot=iplot+1
     endif
     ipict=ipict+1

  endwhile

  for ifile=0,nfile-1 do close,ifile+10
  print
  !p.multi=0
  !p.title=''
  !z.title=''

  if videosave then begin
     videoobject = 0            ; close video file
     print,'Created ',videotime,' sec long video file ',videofile+'.'+savemovie
  endif
  if savemovie eq 'ps' then set_plot,'X'
  ;; Restore velpos array
  velpos=velpos0 & velpos0=0
  if doanimate then xinteranimate,5,/keep_pixmaps
end

;===========================================================================
pro slice_data

;    Animate slices of 3D data in the w array (usually read by 
;    read_data or animate_data).
;    The functions are defined in the "Idl/funcdef.pro" file.
;    "cut" should be defined in 2D, possibly using the grid2d index array.
;
;    Usage:
;
; func='...'     ; functions to plot (optional)
; plotmode='...' ; plot mode(s) (optional)
; slicedir=...   ; select slice direction (1, 2 or 3) (optional)
; slice_data

  common debug_param & on_error, onerror

  common slice_param
  common ask_param
  common getpict_param
  common plot_param
  common plotfunc_param
  common file_head
  common plot_data
  common animate_param

  if not keyword_set(nfile) then read_data

  if nfile gt 1 then begin
     print,'More than one files were read...'
     print,'Probably w is from file ',filenames(nfile-1)
     nfile=1
  endif

  siz=size(x)

  if siz(0) ne 4 then begin
     print,'The x array is not 3 dimensional'
     return
  endif
  n1=siz(1) & n2=siz(2) & n3=siz(3)

  print,'======= CURRENT SLICING PARAMETERS ================'
  print,'firstslice=',firstslice,', dslice=',dslice,$
        ', nslicemax=',nslicemax,', savemovie (y/n)=',savemovie,$
        FORMAT='(a,i4,a,i4,a,i4,a,a)'
  print,'colorlevel=',colorlevel,', contourlevel=',contourlevel,$
        ', velvector=',velvector, $
        FORMAT='(a,i3,a,i3,a,i4)'
  if keyword_set(multiplot) then begin
     siz=size(multiplot)
     if siz(0) eq 0 then begin
        if multiplot gt 0 then multiplot=[multiplot,1,1] $
        else                   multiplot=[1,-multiplot,1]
     endif
     print,'multiplot= ',multiplot,', axistype (coord/cells)=',axistype,$
           ', fixaspect= ',fixaspect,$
           FORMAT='(a,"[",i2,",",i2,",",i2,"]",a,a,a,i1)'
  endif else $
     print,'multiplot= 0 (default), axistype (coord/cells)=',axistype,$
           ', fixaspect= ',fixaspect,$
           FORMAT='(a,a,a,i1)'
  print,'bottomline=',bottomline,', headerline=',headerline,$
        FORMAT='(a,i1,a,i1)'

  if keyword_set(cut) then help,cut
  if keyword_set(velpos) then help,velpos
  velpos0=velpos
  velspeed=0

  print
  help,x
  asknum,'slicedir (1, 2, or 3)',slicedir,doask
  asknum,'firstslice',firstslice,doask
  asknum,'dslice',dslice,doask

  siz=size(x)
  if dslice gt 0 then nslice=(siz(slicedir)-firstslice)/dslice+1 $
  else                nslice=-firstslice/dslice+1
  if nslice gt nslicemax then nslice=nslicemax
  if nslice lt 1 then nslice=0
  print,'Number of slices:',nslice
  if nslice lt 1 then begin
     print,'There are no slices to be shown!'
     print,'   Check firstslice=',firstslice,' and dslice=',dslice
     if siz(slicedir) lt firstslice then $
        print,'   The value of firstslice is larger than the'
     print,'   grid size in slice direction (',slicedir,') =',siz(slicedir)
     retall
  endif

                                ; store 3D data
  var3d = variables
  x3d   = x
  w3d   = w

  case slicedir of
     1:begin
        x=dblarr(n2,n3,2)
        w=dblarr(n2,n3,nw)
        variables=var3d(1:*)
        grid2d=lindgen(n2,n3)
     end
     2:begin
        x=dblarr(n1,n3,2)
        w=dblarr(n1,n3,nw)
        variables=[var3d(0),var3d(2:*)]
        grid2d=lindgen(n1,n3)
     end
     3:begin
        x=dblarr(n1,n2,2)
        w=dblarr(n1,n2,nw)
        variables=[var3d(0:1),var3d(3:*)]
        grid2d=lindgen(n1,n2)
     end
  endcase

  help,grid2d

  print,'======= PLOTTING PARAMETERS ========================='
  read_plot_param

  usereg=0

  if keyword_set(multiplot) then begin
     multix=multiplot(0)
     multiy=multiplot(1)
     nslice1=(multix*multiy)/(nplot*nfile)
     !p.multi=[0,multix,multiy,0,multiplot(2)]
  endif else begin
     multix=long(sqrt(nplot-1)+1)
     multiy=long((nplot-1)/multix+1)
     nslice1=1
     !p.multi=[0,multix,multiy,0,0]
  endelse

  print,'======= DETERMINE PLOTTING RANGES ==================='
  read_limits

  if not noautorange then $
     for islice=1,nslice do begin
     ix=dslice*(islice-1)+firstslice-1
     case slicedir of
        1:begin
           x(*,*,*)=x3d(ix,*,*,1:2)
           w(*,*,*)=w3d(ix,*,*,*)
        end
        2:begin
           x(*,*,0)=x3d(*,ix,*,0)
           x(*,*,1)=x3d(*,ix,*,2)
           w(*,*,*)=w3d(*,ix,*,*)
        end
        3:begin
           x(*,*,*)=x3d(*,*,ix,0:1)
           w(*,*,*)=w3d(*,*,ix,*)
        end
     endcase
     
     first= islice eq 1
     get_limits,first

  endfor

  print
  for ifunc=0,nfunc-1 do $
     print,'Min and max value for ',funcs(ifunc),':',fmin(ifunc),fmax(ifunc)


                                ;==== ANIMATE THE SLICES
  doanimate= nslice gt nslice1 and !d.name eq 'X'
  if !d.name eq 'X' then begin
     if !d.window lt 0 then window
     wshow
  endif
  if doanimate then xinteranimate,set=[!d.x_size,!d.y_size,nslice]

  islice1=0                     ; slice index in a multiplot frame
  iplot=0                       ; plot index for animation
  for islice=1,nslice do begin
     ix=dslice*(islice-1)+firstslice-1
     dirname=var3d(slicedir-1)
     case slicedir of
        1:begin
           x(*,*,*)=x3d(ix,*,*,1:2)
           w(*,*,*)=w3d(ix,*,*,*)
           height=x3d(ix,0,0,0)
           info1='i'+dirname+'='+string(ix,format='(i4)')
        end
        2:begin
           x(*,*,0)=x3d(*,ix,*,0)
           x(*,*,1)=x3d(*,ix,*,2)
           w(*,*,*)=w3d(*,ix,*,*)
           height=x3d(0,ix,0,1)
           info1='i'+dirname+'='+string(ix,format='(i4)')
        end
        3:begin
           x(*,*,*)=x3d(*,*,ix,0:1)
           w(*,*,*)=w3d(*,*,ix,*)
           height=x3d(0,0,ix,2)
           info1='i'+dirname+'='+string(ix,format='(i4)')
        end
     endcase
     info2=dirname+'='+string(height)
     case bottomline of
        0: info=''
        1: info=info2
        2: info=info1+' '+info2
        3: info='time='+string(time,format='(g12.5)')+' '+info1+' '+info2
     endcase

     if not keyword_set(noerase) and islice1 eq 0 then erase

     if velrandom then velpos=0

     rBody3d = rBody
     if abs(height) ge rBody3d then rBody=0.0 $
     else                           rBody=sqrt(rBody3d^2 - height^2)

     plot_func

     xyouts,5+(plotix*!d.x_size)/multix, $
            (1-dyslicelabel)*!d.y_size $
            +  dyslicelabel*(plotiy*!d.y_size)/multiy,/DEV,info
     putheader,1,1,0,0,headerline,headline,nx

     if doanimate then xinteranimate,frame=iplot,window=!d.window
     
     islice1=islice1+1
     if islice1 ge nslice1 then begin
        islice1=0
        iplot=iplot+1
     endif
  endfor

  ;; restore 3d state
  slice_data_restore
  
  print
  !p.multi=0
  !p.title=''

  ;; Restore velpos array
  velpos=velpos0 & velpos0=0

  if doanimate then xinteranimate,5,/keep_pixmaps
end

;=============================================================================
pro slice_data_restore

  ; restore the 3D state after slice_data
  ; this can be useful is slice_data crashed

  common debug_param & on_error, onerror

  common plot_data
  common slice_param
  common file_head
  common plot_param

  x         = x3d
  w         = w3d
  variables = var3d
  rBody     = rBody3d
end

;=============================================================================
pro read_log_data

; Read the log data from 1, 2, or 3 files into the wlog, wlog1, wlog2 arrays
; Store the names of the variables into wlognames, wlognames1 and wlognames2
; Calculate and store log time into the logtime, logtime1 and logtime2
; arrays. The time units are set by timeunit.

  common debug_param & on_error, onerror

  common getlog_param
  common log_data
  common ask_param

  nlogfile=0
  askstr,'logfilename(s) ',logfilename,doask
  string_to_array, logfilename, logfilenames, nlogfile, /wildcard

  if not keyword_set(logfilenames) then begin
     print,'Error in read_log_data: no matching filename was found.'
     retall
  endif

  if nlogfile gt 10 then begin
     print,'Error in read_log_data: cannot handle more than 10 files.'
     print,'nlogfile     = ', nlogfile
     print,'logfilenames = ', logfilenames
     retall
  endif

  get_log,   logfilenames(0), wlog,  wlognames,  logtime, timeunit, /verbose
  if nlogfile ge 2 then $
     get_log, logfilenames(1), wlog1, wlognames1, logtime1, timeunit, verbose='1'
  if nlogfile ge 3 then $
     get_log, logfilenames(2), wlog2, wlognames2, logtime2, timeunit, verbose='2'
  if nlogfile ge 4 then $
     get_log, logfilenames(3), wlog3, wlognames3, logtime3, timeunit, verbose='3'
  if nlogfile ge 5 then $
     get_log, logfilenames(4), wlog4, wlognames4, logtime4, timeunit, verbose='4'
  if nlogfile ge 6 then $
     get_log, logfilenames(5), wlog5, wlognames5, logtime5, timeunit, verbose='5'
  if nlogfile ge 7 then $
     get_log, logfilenames(6), wlog6, wlognames6, logtime6, timeunit, verbose='6'
  if nlogfile ge 8 then $
     get_log, logfilenames(7), wlog7, wlognames7, logtime7, timeunit, verbose='7'
  if nlogfile ge 9 then $
     get_log, logfilenames(8), wlog8, wlognames8, logtime8, timeunit, verbose='8'
  if nlogfile gt 9 then $
     get_log, logfilenames(9), wlog9, wlognames9, logtime9, timeunit, verbose='9'

end

;=============================================================================
pro plot_log_data

  common debug_param & on_error, onerror

  common ask_param
  common plotlog_param
  common log_data

  if not keyword_set(wlog) then read_log_data

  askstr,'logfunc(s)     ', logfunc, doask

  if !d.name eq 'X' then begin
     if !d.window lt 0 then window
     wshow
  endif

  plot_log

end

;=============================================================================
function reform2,x

  common debug_param & on_error, onerror

  ;; Remove all degenerate dimensions from x

  if n_elements(x) lt 2 then return,x

  siz = size(x)
  siz = siz(1:siz(0))
  return, reform(x, siz(where(siz gt 1)))

end

;=============================================================================
function log_time,wlog,wlognames,timeunit

; Obtain time in hours from wlog and wlognames
; If the log file contains 't' or 'time', simply convert seconds to hours
; If the log file contains date and time, calculate time relative to
; the beginning of the first day. 
; This algorithm works only if the whole file is in the same month.

  common debug_param & on_error, onerror

  nwlog = n_elements(wlognames)
  hours = 0.0*wlog(*,0)

  if nwlog lt 0 then return, hours

  itime = -1
  istep = -1
  iyear = -1
  iday  = -1
  ihour = -1
  imin  = -1
  isec  = -1
  imsc  = -1
  for i = 0, nwlog-1 do begin
     varname = strlowcase(wlognames(i))
     case varname of
        'time'       : itime = i
        't'          : if wlognames(i) eq 't' then itime=i
        'step'       : istep = i
        'it'         : istep = i
        'year'       : iyear = i
        'yr'         : iyear = i
        'yy'         : iyear = i
        'day'        : iday  = i
        'dy'         : iday  = i
        'dd'         : iday  = i
        'hour'       : ihour = i
        'hr'         : ihour = i
        'hh'         : ihour = i
        'minute'     : imin  = i
        'min'        : imin  = i
        'mn'         : imin  = i
        'mm'         : imin  = i
        'second'     : isec  = i
        'sec'        : isec  = i
        'sc'         : isec  = i
        'ss'         : isec  = i
        'millisecond': imsc  = i
        'millisec'   : imsc  = i
        'msecond'    : imsc  = i
        'msec'       : imsc  = i
        'msc'        : imsc  = i
        else:
     endcase
  endfor

  if itime gt -1 then begin
     hours = wlog(*,itime)/3600.0
  endif else begin
     if iyear eq -1 or iday eq -1 then begin
        if ihour gt -1 then begin
           hours = wlog(*,ihour)
           for i=1L, n_elements(hours)-1 do $
              while hours(i) lt hours(i-1) do $
                 hours(i) = hours(i) + 24.0
        endif
     endif else begin
        ; calculate time from year, day, hour, minute, second, millisec
        nday = 0
        daylast  = wlog(0,iday)
        for i = 0L, n_elements(hours) - 1 do begin
           dday = wlog(i,iday) - daylast
           if dday ne 0 then begin
              ; if we go to the next month we assume
              ; we ended on the last day of the month(?)
              if dday lt 0 then dday = wlog(i,iday)
              daylast = wlog(i,iday)
              nday = nday + dday
           endif
           hours[i] = nday*24.0
           if ihour eq iday + 1 then hours[i] = hours[i] + wlog(i,ihour)
           if imin  eq iday + 2 then hours[i] = hours[i] + wlog(i,imin)/60.0
           if isec  eq iday + 3 then hours[i] = hours[i] + wlog(i,isec)/3600.0
           if imsc  eq iday + 4 then hours[i] = hours[i] + wlog(i,imsc)/3.6e6
        endfor
     endelse
  endelse

  if max(hours) eq min(hours) then begin
     if istep gt -1 then begin
        print, 'Could not find time information, using steps'
        hours = wlog(*,istep)
     endif else begin
        print, 'Could not find time information, using line number'
        hours = findgen(n_elements(wlog(*,0)))
     endelse
  endif

  logtime = hours

  if n_elements(timeunit) gt 0 then begin
     case timeunit of
        'd': logtime = hours/24
        '1': logtime = hours*3600
        's': logtime = hours*3600
        'm': logtime = hours*60
        'millisec': logtime = hours*3600e3
        'microsec': logtime = hours*3600e6
        'ns'    : logtime = hours*3600e9
        else: 
     endcase
  endif

  return, logtime
end

;=============================================================================
function log_func, wlog, varnames, varname, error

  common debug_param & on_error, onerror

  error = 0
  ivar  = where(strlowcase(varnames) eq strlowcase(varname)) & ivar = ivar(0)
                                ; Variable is found, return with array
  if ivar ge 0 then return, wlog(*,ivar)

                                ; Try calculating temperature or pressure
  if varname eq 'T' then begin
                                ; Convert p[nPa]/n[/cc] to T[eV]
     ivar = where( varnames eq 'p')   & ivar = ivar(0)
     jvar = where( varnames eq 'rho') & jvar = jvar(0)
     if ivar ge 0 and jvar ge 0 then $
        return,6241.5*wlog(*,ivar)/wlog(*,jvar)

  endif else if varname eq 'p' then begin
                                ; Convert T[eV]*n[/cc] to p[nPa] to T[eV]
     ivar = where( varnames eq 'T')   & ivar = ivar(0)
     jvar = where( varnames eq 'rho') & jvar = jvar(0)
     if ivar ge 0 and jvar ge 0 then $
        return, wlog(*,ivar)*wlog(*,jvar)/6241.5
  endif

  error = 1
  return,0*wlog(*,0)
  
end

;=============================================================================
pro black_background

  common debug_param & on_error, onerror

  !p.background =   0
  !p.color      = 255
  POLYFILL, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=0 ; draw black box

end
;=============================================================================
pro white_background

  common debug_param & on_error, onerror

  !p.background = 255
  !p.color      =   0

end
;=============================================================================
pro open_file,unit,filename,filetype

   common debug_param & on_error, onerror

   i = strpos(filename,'.gz')
   if i gt 0 then begin
       filenamenew = strmid(filename,0,i)
       print,'gunzip -c '+filename+' > '+filenamenew
       spawn,'gunzip -c '+filename+' > '+filenamenew
       filename = filenamenew
   endif

   close,unit
   case filetype of
       'log'   :openr,unit,filename
       'ascii' :openr,unit,filename
       'binary':openr,unit,filename,/f77_unf
       'real4' :openr,unit,filename,/f77_unf
       'BINARY':openr,unit,filename,/f77_unf
       'REAL4' :openr,unit,filename,/f77_unf
       'IPIC3D':
       else    :print,'open_file: unknown filetype:',filetype
   endcase
end

;=============================================================================
pro get_file_types

  common debug_param & on_error, onerror

  common getpict_param

  filetypes    = strarr(nfile)
  npictinfiles = intarr(nfile)
  for ifile=0, nfile-1 do begin
     l = strlen(filenames(ifile)) - 4
     if   strpos(filenames(ifile),'.log') eq l $
        or strpos(filenames(ifile),'.sat') eq l then begin
        filetypes(ifile)    = 'log'
        npictinfiles(ifile) = 1
     endif else begin
        if strpos(filenames(ifile),'setting') ge 0 then begin
                                ; For example if 
                                ; filenames(ifile) = 'output/settings_region0.hdf'   then 
                                ; dirname          = 'output/'                       and
                                ; regionname       =                '_region0.hdf'.                      
           ibegin     = strpos(filenames(ifile), 'settings')
           dirname    = strmid(filenames(ifile), 0, ibegin)

           regionname ='.hdf'
           if strpos(filenames(ifile),'settings_region') ge 0 then begin
              regionname = strmid(filenames(ifile), ibegin  + strlen('settings'))
           endif
           file_id    = H5F_OPEN(dirname+'proc0'+regionname)
           group_id   = H5G_OPEN(file_id, '/fields/Bx')
           npictinfiles(ifile) = H5G_GET_NUM_OBJS(group_id)
           h5G_CLOSE, group_id
           H5F_CLOSE, file_id
           filetypes(ifile)   ='IPIC3D'

        endif else begin
           ;; Obtain filetype based on the length info in the first 4 bytes
           close,10
           openr,10,filenames(ifile)
           lenhead=long(1)
           readu,10,lenhead
           if lenhead ne 79 and lenhead ne 500 then ftype='ascii' else begin
              ;; The length of the 2nd line decides between real4 and binary
              ;; since it contains the time, which is real*8 or real*4
              head=bytarr(lenhead+4)
              len=long(1)
              readu,10,head,len
              case len of
                 20: ftype='real4'
                 24: ftype='binary'
                 else: begin
                    print,'Error in get_file_types: strange unformatted file:',$
                          filenames(ifile)
                    retall
                 end
              endcase
              if lenhead eq 500 then ftype = strupcase(ftype)
           endelse
           close,10
           
           ;; Obtain file size and number of snapshots
           open_file,1,filenames(ifile),ftype
           status=fstat(1)
           fsize=status.size

           pointer=long64(0)
           pictsize=long64(1)
           ipict=0
           while pointer lt fsize do begin
                                ; Obtain size of a single snapshot
              point_lun,1,pointer
              get_file_head,1,filenames(ifile),ftype,pictsize=pictsize
              ipict   = ipict+1
              pointer = pointer + pictsize
           endwhile
           close,1

           npictinfiles(ifile)=ipict
           filetypes(ifile)   =ftype
        endelse
     endelse
  endfor
end

;=============================================================================
pro show_head, ifile

  common debug_param & on_error, onerror

  common ask_param
  common getpict_param
  common file_head

  print,   'filename   = ',filenames(ifile), format="(a,a)"
  print,   'filetype   = ',filetypes(ifile), format="(a,a)"
  print,   'headline   = ',strtrim(headline,2), format="(a,a)"
  print,   'it         = ',it,       format="(a,i8)"
  print,   'time       = ',time,     format="(a,g15.8)"
  print,   'gencoord   = ',gencoord, format="(a,i8)"
  print,   'ndim       = ',ndim,     format="(a,i8)"
  print,   'neqpar     = ',neqpar,   format="(a,i8)"
  print,   'nw         = ',nw,       format="(a,i8)"
  print,   'nx         = ',nx,       format="(a,3i8)"
  if neqpar gt 0 then $
     print,'parameters = ',eqpar,    format="(a,100g15.8)"
  print,   'coord names=',variables(0:ndim-1)
  print,   'var   names=',variables(ndim:ndim+nw-1)
  if neqpar gt 0 then $
     print,'param names=',variables(ndim+nw:*)

end
;=============================================================================
pro get_file_head, unit, filename, filetype, pictsize=pictsize

  common debug_param & on_error, onerror

  common file_head

  ftype = strlowcase(filetype)

  if ftype eq filetype then lenstr = 79 else lenstr = 500

  ;; Type definitions
  headline=''
  for i=1, lenstr do headline=headline+' '
  it=long(1)
  ndim=long(1)
  neqpar=long(0)
  eqpar=0.0
  nw=long(1)
  varname=''
  for i=1, lenstr do varname=varname+' '

  ;; Remember pointer position at beginning of header
  if ftype ne 'ipic3d' then point_lun,-unit,pointer0

  ;; Read header
  case ftype of
     'ipic3d':begin
        tmppict = 0
        tmperror = 0
        get_pict_hdf, filename, tmppict, tmperror, 0
     end
     'log': begin
        readf,unit,headline
        readf,unit,varname
        nw=n_elements(strsplit(varname))-2
        varname ='hour '+varname
                                ; reset pointer
        point_lun,unit,pointer0
        it=0
        time=0.0
        gencoord=0
        ndim=1
        nx=lonarr(1)
        nx(0)=1
     end
     'ascii': begin
        time=double(1)
        readf,unit,headline
        readf,unit,it,time,ndim,neqpar,nw
        gencoord=(ndim lt 0)
        ndim=abs(ndim)
        nx=lonarr(ndim)
        readf,unit,nx
        if neqpar gt 0 then begin
           eqpar=dblarr(neqpar)
           readf,unit,eqpar
        endif
        readf,unit,varname
     end
     'binary':begin
        time=double(1)
        readu,unit,headline
        readu,unit,it,time,ndim,neqpar,nw
        gencoord=(ndim lt 0)
        ndim=abs(ndim)
        nx=lonarr(ndim)
        readu,unit,nx
        if neqpar gt 0 then begin
           eqpar=dblarr(neqpar)
           readu,unit,eqpar
        endif
        readu,unit,varname
     end
     'real4': begin
        time=float(1)
        readu,unit,headline
        readu,unit,it,time,ndim,neqpar,nw
        gencoord=(ndim lt 0)
        ndim=abs(ndim)
        nx=lonarr(ndim)
        readu,unit,nx
        if neqpar gt 0 then begin
           eqpar=fltarr(neqpar)
           readu,unit,eqpar
        endif
        readu,unit,varname
     end
     else: begin
        print,'get_file_head: unknown filetype',filetype
        retall
     end
  endcase

  if keyword_set(pictsize) then begin
                                ; Calculate the picture size
                                ; Header length
     point_lun,-unit,pointer1
     headlen=pointer1-pointer0
                                ; Number of cells
     nxs=long64(1)
     for idim=1,ndim do nxs=nxs*nx(idim-1)
                                ; Snapshot size = header + data + recordmarks
     case ftype of
        'log'   :pictsize = 1
        'ascii' :pictsize = headlen + (18*(ndim+nw)+1)*nxs
        'binary':pictsize = headlen + 8*(1+nw)+8*(ndim+nw)*nxs
        'real4' :pictsize = headlen + 8*(1+nw)+4*(ndim+nw)*nxs
     endcase
  endif

  ;; Set variables array
  string_to_array,varname,variables,nvar,/arraysyntax
end

;=============================================================================
pro get_pict_hdf, filename, npict, error, getdata

  common debug_param & on_error, onerror

  common file_head
  common plot_data

    ;;;;;;;;;;;;;;;;; SIM PARAMETERS ;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ibegin     = strpos(filename, 'settings')
  dirname    = strmid(filename, 0, ibegin)
  regionname ='.hdf'
  if strpos(filename,'settings_region') ge 0 then begin
     regionname = strmid(filename, ibegin  + strlen('settings'))
  endif

  Param = H5_PARSE(filename) 

  nxyz_D = [Param.COLLECTIVE.NXC._DATA(0),Param.COLLECTIVE.NYC._DATA(0),$
            Param.COLLECTIVE.NZC._DATA(0)]    

  gencoord = 0
  headline = 'PIC CGS units'

  idims = where(nxyz_D GT 1, ndim)
  nx = nxyz_D(idims)
  ;; Bx By Bz Ex Ey Ez + [rho] + [Jx,Jy,Jz]. rho, Jx, Jy, Jz are optional.
  nw =10 + Param.COLLECTIVE.NS._DATA(0)*10

  eqpar = [Param.COLLECTIVE.BX0._DATA(0),Param.COLLECTIVE.BY0._DATA(0),$
           Param.COLLECTIVE.BZ0._DATA(0)]
  neqpar= 3
  
  dxyz_D = [Param.COLLECTIVE.Dx._DATA(0),$
            Param.COLLECTIVE.Dy._DATA(0),$
            Param.COLLECTIVE.Dz._DATA(0)]   

  if getdata then begin
     case ndim of
        1:begin
           x=DBLARR(nx(0),ndim)
           w=DBLARR(nx(0),nw)
           for x0=0L,nx(0)-1 do begin
              x(x0,0:ndim-1) = x0*dxyz_D(0:nDim-1)
           endfor
        end
        2:begin
           x=DBLARR(nx(0),nx(1),ndim)
           w=DBLARR(nx(0),nx(1),nw)
           for x1=0L,nx(1)-1 do begin
              for x0=0L,nx(0)-1 do begin
                 x(x0,x1,0:ndim-1) = [x0,x1]*dxyz_D(0:nDim-1)
              endfor
           endfor
        end
        3:begin
           x=DBLARR(nx(0),nx(1),nx(2),ndim)
           w=DBLARR(nx(0),nx(1),nx(2),nw)
           for x2=0L,nx(2)-1 do begin
              for x1=0L,nx(1)-1 do begin
                 for x0=0L,nx(0)-1 do begin
                    x(x0,x1,x2,0:ndim-1) = [x0,x1,x2]*dxyz_D(0:nDim-1)
                 endfor
              endfor
           endfor
        end
     endcase
  endif

  ;; processor layout
  nproc = Param.TOPOLOGY.Nprocs._DATA(0)

  Proc_D = [Param.TOPOLOGY.XLEN._DATA(0),$
            Param.TOPOLOGY.YLEN._DATA(0),$
            Param.TOPOLOGY.ZLEN._DATA(0)]

  nExtra_D     = INTARR(3)
  iCrowd_D     = INTARR(3)
  nxyzLocal_PD = INTARR(nproc,3)
  nxyzLocal0_D = nxyz_D/Proc_D
  nCrowded_D   = nxyz_D - nxyzLocal0_D*Proc_D
  
  ;; index range
  MinIJK_PD = INTARR(nproc,3)
  MaxIJK_PD = INTARR(nproc,3)
  iproc=0
  for ip=0,Proc_D(0)-1 do begin
     for jp=0,Proc_D(1)-1 do begin
        for kp=0,Proc_D(2)-1 do begin
           iCrowd_D(*) = 0
           nExtra_D = nCrowded_D
           if ip LT nCrowded_D(0) then begin
              iCrowd_D(0) = 1
              nExtra_D(0) = ip
           endif
           if jp LT nCrowded_D(1) then begin
              iCrowd_D(1) = 1
              nExtra_D(1) = jp
           endif
           if kp LT nCrowded_D(2) then begin
              iCrowd_D(2) = 1
              nExtra_D(2) = kp 
           endif
           nxyzLocal_PD(iproc,*) = nxyzLocal0_D + iCrowd_D           
           MinIJK_PD(iproc,*) = [ip,jp,kp]*nxyzLocal0_D + nExtra_D
           MaxIJK_PD(iproc,*) = MinIJK_PD(iproc,*) + nxyzLocal_PD(iproc,*) - 1
           iproc=iproc+1
        endfor
     endfor
  endfor

  MinIJK_PD = MinIJK_PD(0:nproc-1,idims)
  MaxIJK_PD = MaxIJK_PD(0:nproc-1,idims)

  ;;;;;;;;;;;;;;;;; GETTING TIMELINE ++ ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  filename= dirname+'proc0' + regionname
  file_id = H5F_OPEN(filename)
  
  ;; Find the cronological timeline index SortIdx_I
  group_id = H5G_OPEN(file_id, '/fields/Bx')
  nObj = H5G_GET_NUM_OBJS(group_id)
  Step_I = LONARR(nObj)         ;; store time index of saved snapshots
  for iObj=0,nObj-1 do begin
     ObjName = H5G_GET_OBJ_NAME_BY_IDX(group_id,iObj)
     Step_I(iObj) = STRMID(ObjName,6)
  endfor
  h5G_CLOSE, group_id
  H5F_CLOSE, file_id
  SortIdx_I = SORT(Step_I)
  
  ;;bounding npict
  if npict lt 1 then npict=1
  if npict gt nObj then begin
     error = 1
     return
  endif

  ;; Find iteration and time 
  it= Step_I[SortIdx_I(npict-1)]
  time =  it*Param.COLLECTIVE.Dt._DATA(0)
  
  ;; Setting up "variables"
  nVar = nw +ndim + neqpar
  variables = STRARR(nVAr)
  DimName = ['x','y','z']
  variables(0:ndim-1) = DimName(0:ndim-1)
  variables(nVar-neqpar:nVar-1) = ['B0x','B0y','B0z']
  
  ;;;;;;;;;;;;;;;;; DATA GATHERING ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; loop over all files
  for iproc=0,nproc-1 do begin
     
     filename = dirname+'proc' + string(iproc,FORMAT='(I0)') + regionname

     file_id = H5F_OPEN(filename)
     iMin = MinIJK_PD(iproc,0)
     iMax = MaxIJK_PD(iproc,0)
     nxyz_D  = nxyzLocal_PD(iproc,*)

     if ndim gt 1 then begin 
        jMin = MinIJK_PD(iproc,1)
        jMax = MaxIJK_PD(iproc,1)
     end

     if ndim gt 2 then begin
        kMin = MinIJK_PD(iproc,2)
        kMax = MaxIJK_PD(iproc,2)
     end

     iw = 0
     ; Go through all field variables
     group_id = H5G_OPEN(file_id, '/fields')
     nFields = H5G_GET_NUM_OBJS(group_id)
     for iFields=0, nFields-1 do begin
        get_hdf_pict,group_id,iFields, SortIdx_I(npict-1), nxyz_D,$
                     -1, pict, varname, getdata

        if getdata then begin
           case ndim of
              1:w(iMin:iMax,iw) = pict
              2:w(iMin:iMax,jMin:jMax,iw) = pict
              3:w(iMin:iMax,jMin:jMax,kMin:kMax,iw) = pict
           endcase
        end
        variables(ndim+iw) = varname
        iw = iw +1
     endfor
     h5G_CLOSE, group_id  

     ;; Get all moments of the distribution function
     group_id = H5G_OPEN(file_id, '/moments')
     nVarable = H5G_GET_NUM_OBJS(group_id)     

     iSpecies = 0
     for iVar=0,nVarable-1 do begin
                                ; Read rho Jx Jy Jz...
        VarName = H5G_GET_OBJ_NAME_BY_IDX(group_id,iVar)
        if strmid(VarName,0,7) ne 'species' then begin
           get_hdf_pict,group_id,iVar,SortIdx_I(npict-1),nxyz_D,$
                        -1,pict,varname,getdata           
           if getdata then begin
              case ndim of
                 1:w(iMin:iMax,iw) = pict
                 2:w(iMin:iMax,jMin:jMax,iw) = pict
                 3:w(iMin:iMax,jMin:jMax,kMin:kMax,iw) = pict
              endcase
           end
           variables(ndim+iw) = varname
           iw = iw +1
        endif else begin
           SpeciesName = H5G_GET_OBJ_NAME_BY_IDX(group_id,iVar)
           Moment_id = H5G_OPEN(group_id, SpeciesName)
           nMoment = H5G_GET_NUM_OBJS(Moment_id)
           for iMoment=0, nMoment-1 do begin
              get_hdf_pict,Moment_id,iMoment,SortIdx_I(npict-1),nxyz_D,$
                           iSpecies,pict,varname,getdata
              if getdata then begin
                 case ndim of
                    1:w(iMin:iMax,iw) = pict
                    2:w(iMin:iMax,jMin:jMax,iw) = pict
                    3:w(iMin:iMax,jMin:jMax,kMin:kMax,iw) = pict
                 endcase
              end
              variables(ndim+iw) = varname
              iw = iw +1
           endfor
           iSpecies = iSpecies + 1
           h5G_CLOSE, Moment_id  
        endelse 
     endfor
     h5G_CLOSE, group_id  
     H5F_CLOSE, file_id
  endfor
end 

;=============================================================================
pro get_hdf_pict,group_id,iGroup,ipict,nx,iSpecies,pictout,name,getdata

  common debug_param & on_error, onerror

  GroupName = H5G_GET_OBJ_NAME_BY_IDX(group_id,iGroup)
  name= GroupName
  if iSpecies ge 0 then $
     name= GroupName +"S"+STRING(iSpecies,FORMAT='(I02)')
  moment_id = H5G_OPEN(group_id, GroupName)
  ; get field data form npict time pict
  pictname = H5G_GET_OBJ_NAME_BY_IDX(moment_id,ipict)
  if getdata then begin
     pict_id = H5D_OPEN(moment_id,pictname)
     pict = H5D_READ(pict_id)
     pict = 0.5*(pict(0:nx(2)-1,*,*) + pict(1:nx(2),*,*))
     pict = 0.5*(pict(0:nx(2)-1,0:nx(1)-1,*) + pict(0:nx(2)-1,1:nx(1),*))
     pict = 0.5*(pict(0:nx(2)-1,0:nx(1)-1,0:nx(0)-1) + pict(0:nx(2)-1,0:nx(1)-1,1:nx(0)))
     pictout = reform(TRANSPOSE(pict(0:nx(2)-1,0:nx(1)-1,0:nx(0)-1),[2,1,0]))
     H5D_CLOSE, pict_id
  endif
  h5G_CLOSE, moment_id

end
;=============================================================================
pro get_pict, unit, filename, filetype, npict, error
  
  common debug_param & on_error, onerror

  if filetype eq 'IPIC3D' then begin 
     error=0
     get_pict_hdf, filename, npict, error, 1

  endif else begin

     error=0

     if(eof(unit))then begin
        error=1
        return
     endif

                                ; Get current pointer position
     point_lun,-unit,pointer

                                ; Skip npict-1 snapshots
     ipict=0
     pictsize=1
     while ipict lt npict-1 and not eof(unit) do begin
        ipict=ipict+1
        get_file_head,unit,filename,filetype,pictsize=pictsize
        pointer=long64(pointer) + pictsize
        point_lun,unit,pointer
     endwhile

                                ; Backup 1 snapshot if end of file
     if eof(unit) then begin
        error=1
        point_lun,unit,pointer-pictsize
     endif

                                ; Read header information
     get_file_head, unit, filename, filetype

                                ; Read data
     case strlowcase(filetype) of
        'log':    get_pict_log ,unit
        'ascii':  get_pict_asc ,unit, npict
        'binary': get_pict_bin ,unit, npict
        'real4':  get_pict_real,unit, npict
        else:    begin
           print,'get_pict: unknown filetype:',filetype
           error=1
           close,unit
        end
     endcase

  endelse

  set_units

end

;=============================================================================
pro get_pict_log, unit

  common debug_param & on_error, onerror

  common plot_data
  common file_head

  get_log, unit, w, wlognames, x, 's'

  ndim = 1
  nx(0)= n_elements(x)
  nw   = n_elements(wlognames)
  
end

;=============================================================================
pro get_pict_asc,unit,npict

  common debug_param & on_error, onerror

  common plot_data
  common file_head

  ;----------------------------------------
  ; Read coordinates and values row by row
  ;----------------------------------------
  wrow=dblarr(nw)
  xrow=dblarr(ndim)
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    x=dblarr(nx(0),ndim)
    w=dblarr(nx(0),nw)
    for x0=0L,nx(0)-1 do begin
      readf,unit,xrow,wrow
      x(x0,0:ndim-1)=xrow(0:ndim-1)
      w(x0,0:nw-1)  =wrow(0:nw-1)
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    x=dblarr(nx(0),nx(1),ndim)
    w=dblarr(nx(0),nx(1),nw)
    for x1=0L,nx(1)-1 do begin
      for x0=0L,nx(0)-1 do begin
        readf,unit,xrow,wrow
        x(x0,x1,0:ndim-1)=xrow(0:ndim-1)
        w(x0,x1,0:nw-1)  =wrow(0:nw-1)
      endfor
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    x=dblarr(nx(0),nx(1),nx(2),ndim)
    w=dblarr(nx(0),nx(1),nx(2),nw)
    for x2=0L,nx(2)-1 do begin
      for x1=0L,nx(1)-1 do begin
        for x0=0L,nx(0)-1 do begin
          readf,unit,xrow,wrow
          x(x0,x1,x2,0:ndim-1)=xrow(0:ndim-1)
          w(x0,x1,x2,0:nw-1)=wrow(0:nw-1)
        endfor
      endfor
    endfor
  end
  endcase
end

;=============================================================================
pro get_pict_bin,unit,npict

  common debug_param & on_error, onerror

  common plot_data
  common file_head

  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=dblarr(n1,ndim)
    w=dblarr(n1,nw)
    wi=dblarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=dblarr(n1,n2,ndim)
    w=dblarr(n1,n2,nw)
    wi=dblarr(n1,n2)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=dblarr(n1,n2,n3,ndim)
    w=dblarr(n1,n2,n3,nw)
    wi=dblarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;=============================================================================
pro get_pict_real,unit,npict

  common debug_param & on_error, onerror

  common file_head
  common plot_data

  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=fltarr(n1,ndim)
    w=fltarr(n1,nw)
    wi=fltarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=fltarr(n1,n2,ndim)
    w=fltarr(n1,n2,nw)
    readu,unit,x
    wi=fltarr(n1,n2)
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=fltarr(n1,n2,n3,ndim)
    w=fltarr(n1,n2,n3,nw)
    wi=fltarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;=============================================================================
pro asknum,prompt,var,doask

   common debug_param & on_error, onerror

   if var eq 0 then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+strtrim(string(var),2)+' ? ',tmp
         if tmp ne '' then reads,tmp,var
      endif else print,prompt,'=',var
   endelse
end

;=============================================================================
pro askstr,prompt,var,doask

   common debug_param & on_error, onerror

   if var eq '' then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+var+' ? ',tmp
         if tmp ne '' then var=tmp
      endif else print,prompt,'=',var
   endelse
end

;=============================================================================
pro string_to_array,s,a,n,sep,arraysyntax=arraysyntax, wildcard=wildcard

; If s is an array copy it to a.
; If s is a string then split string s at the sep characters 
;      (default is space) into array a. 
;      Array syntax is expanded: "Var(6:10:2)" --> ["Var06","Var08","Var10"] 
; If n is 0, it will be set to the size of the array on output
; If n is defined, fill up the rest of the array with the last defined element
; If n is defined but smaller than the number of elements in s, print a warning

common debug_param & on_error, onerror

if keyword_set(wildcard) and stregex(s, '[?*[]', /boolean) then begin
   spawn,'/bin/ls '+s, a
   n = n_elements(a)
   return
endif

if n_elements(s) gt 1 then     $
    a0 = s                     $
else if keyword_set(sep) then  $
   a0=strsplit(s,sep,/extract) $
else                           $
   a0=strsplit(s,/extract)

n0=n_elements(a0)

if keyword_set(arraysyntax) then begin
    a1 = strarr(300)
    n1 = 0
    for i0 = 0, n0-1 do begin
 
        s1 = a0(i0)

        ; copy the element into a1 in case it is not array-syntax
        a1(n1) = s1
        n1 = n1 + 1
                                ; check for closing paren 
        if strpos(s1,')') eq strlen(s1)-1 then begin
                                ; chop it off
            s1 = strmid(s1,0,strlen(s1)-1)
                                ; split at opening paren
            s1 = strsplit(s1,'(',/extract)
            if n_elements(s1) eq 2 then begin
                s2 = s1(1)
                s1 = s1(0)
                                ; split second part at colons
                s2 = strsplit(s2,':',/extract)
                                ; set defaults for minimum and stride
                imin   = 1
                stride = 1
                                ; set index for maximum value (0 or 1)
                npart  = n_elements(s2)
                iimax = npart-1 < 1
                                ; set maximum value
                imax  = fix(s2(iimax))
                                ; set minimum value if present
                if npart gt 1 then imin   = fix(s2(0))
                                ; set stride if present
                if npart gt 2 then stride = fix(s2(2))

                                ; create format string
                width   = string(strlen(s2(iimax)),format='(i1)')
                formstr = "(a,i"+width+"."+width+")"

                                ; add new elements to a1 array
                n1 = n1 - 1
                for i = imin,imax,stride do begin
                    a1(n1) = string(s1,i,format=formstr)
                    n1 = n1 + 1
                endfor
            endif
        endif
    endfor
    a0 = a1(0:n1-1)
    n0 = n1
endif

if not keyword_set(n) then begin
   a=a0
   n=n0
endif else if n ge n0 then begin
   a=strarr(n)
   a(0:n0-1)=a0
   if n0 lt n then a(n0:n-1)=a0(n0-1)
endif else begin
   a=strarr(n)
   a=a0(0:n-1)
   print,'Warning: more than',n,' values defined by string: ',s,$
       FORMAT='(a,i3,a,a)'
endelse

end
;=============================================================================
pro arr2arr,a,n

; If a is a scalar or has fewer elements than n then 
; extend it to an array of n elements.
; If a has more elements than n, reduce the number of elements to n.

common debug_param & on_error, onerror

k = n_elements(a)
if      k gt n then a = a(0:n-1) $                    ; truncate array
else if k lt n then for i = k, n-1 do a = [a, a(k-1)] ; extend array

end
;===========================================================================
pro read_plot_param

  common ask_param, doask
  common plot_param
  common plotfunc_param
  common file_head, ndim

  common debug_param & on_error, onerror

  ;; Determine dimension of plots based on cut or ndim,
  ;; calculate reduced cut0 array by eliminating degenerate dimensions
  if keyword_set(cut) then begin
     cut0=reform2(cut)
     siz=size(cut0)
     plotdim=siz(0)
  endif else begin
     plotdim=ndim
     cut0=0
  endelse

  askstr,'func(s) (e.g. rho p ux;uz bx+by -T) ',func,doask
  if plotdim eq 1 then begin
     print,'1D plotmode: plot/plot_io/plot_oi/plot_oo'
     print,'1D +options: max,mean,log,noaxis,over,#c999,#ct999'
     askstr,'plotmode(s)                ',plotmode,doask
     if strmid(plotmode,0,4) ne 'plot' then plotmode='plot'
  endif else begin
     if strmid(plotmode,0,4) eq 'plot' then plotmode='contbar'
     print,'2D plotmode: shade/surface/cont/tv/polar(rad|deg|hour)/velovect/vector/stream'
     print,'2D +options: bar,body,fill,grid,irr,label,max,mean,log'
     print,'2D +options: map,mesh,noaxis,over,usa,white,#c999,#ct999'
     askstr,'plotmode(s)                ',plotmode,doask
  endelse
  askstr,'plottitle(s) (e.g. B [G];J)',plottitle,doask
  askstr,'autorange(s) (y/n)         ',autorange,doask

  nfunc=0
  string_to_array,func,funcs,nfunc
  string_to_array,plotmode,plotmodes,nfunc
  string_to_array,plottitle,plottitles,nfunc,';'
  string_to_array,autorange,autoranges,nfunc

  nplot = nfunc
  funcs1=strarr(nfunc)
  funcs2=strarr(nfunc)
  for ifunc=0,nfunc-1 do begin
     func12=strsplit(funcs(ifunc),';',/extract)
     funcs1(ifunc)=func12(0)
     if n_elements(func12) eq 2 then funcs2(ifunc)=func12(1)
     if strpos(plotmodes(ifunc),'over') ge 0 then nplot=nplot-1
  endfor

end
;===========================================================================
pro read_transform_param

  common debug_param & on_error, onerror

  common ask_param
  common file_head
  common transform_param
  common vector_param
  common plot_data

  if (gencoord or transform eq 'unpolar') and ndim eq 2 then begin
      if transform eq '' then begin
        transform='none'
        askstr,"transform (r=regular/p=polar/u=unpolar/m=my/n=none)",$
          transform,1
      endif else $
        askstr,"transform (r=regular/p=polar/u=unpolar/m=my/n=none)",$
        transform,doask

      ; Complete name
      case transform of
          'r': transform='regular'
          'p': transform='polar'
          'u': transform='unpolar'
          'm': transform='my'
          'n': transform='none'
         else:
      endcase
      ; Get transformation parameters and calculate grid
      case 1 of
        transform eq 'regular':begin
           print,'Generalized coordinates, dimensions for regular grid'
           if n_elements(nxreg) ne 2 then nxreg=[0,0]
           if n_elements(xreglimits) ne 4 then xreglimits=dblarr(4) $
           else xreglimits=double(xreglimits)
           nxreg0=nxreg(0)
           nxreg1=nxreg(1)
           asknum,'nxreg(0) (use negative sign to limit x)',nxreg0,doask
           if nxreg0 lt 0 then begin
               nxreg0=abs(nxreg0)
               xmin=0 & xmax=0
               asknum,'xreglimits(0) (xmin)',xmin,doask
               asknum,'xreglimits(2) (xmax)',xmax,doask
               xreglimits(0)=xmin
               xreglimits(2)=xmax
           endif
           asknum,'nxreg(1) (use negative sign to limit y)',nxreg1,doask
           if nxreg1 lt 0 then begin
               nxreg1=abs(nxreg1)
               ymin=0 & ymax=0
               asknum,'xreglimits(1) (ymin)',ymin,doask
               asknum,'xreglimits(3) (ymax)',ymax,doask
               xreglimits(1)=ymin
               xreglimits(3)=ymax
           endif
           grid=lindgen(nxreg0,nxreg1)
           nxreg=[nxreg0,nxreg1]
           wregpad=0
        end
        transform eq 'polar' or transform eq 'unpolar':begin
            asknum,'Number of vector variables',nvector,doask
            getvectors,nvector,vectors
            grid=lindgen(nx(0),nx(1))
        end
        transform eq 'none':grid=lindgen(nx(0),nx(1))
        else: print,'Unknown value for transform:',transform
      endcase
   endif else if (gencoord or transform eq 'unpolar') and ndim eq 3 then begin
       if transform eq '' then begin
           transform="none"
           askstr,"transform (p=polar/s=sphere/u=unpolar/m=my/n=none)",$
             transform,1
       endif else $
         askstr,"transform (p=polar/s=sphere/u=unpolar/m=my/n=none)",$
         transform,doask
       case transform of
           'p': transform='polar'
           's': transform='sphere'
           'u': transform='unpolar'
           'm': transform='my'
           'n': transform='none'
           else:
       endcase
       case 1 of
           transform eq 'polar' or transform eq 'unpolar' or $
             transform eq 'sphere' :begin
               getvectors,nvector,vectors
               grid=lindgen(nx(0),nx(1),nx(2))
           end
           transform eq 'none': grid=lindgen(nx(0),nx(1),nx(2))
           else: print,'Unknown value for transform:',transform
       endcase
   endif else case ndim of
       1: grid=lindgen(nx(0))
       2: grid=lindgen(nx(0),nx(1))
       3: grid=lindgen(nx(0),nx(1),nx(2))
   endcase

   ;===== GRID HELPS TO CREATE A CUT, E.G.: cut=grid(*,4)

   help,grid
end

;===========================================================================
pro getvectors,nvector,vectors

  common debug_param & on_error, onerror


  if nvector eq 0 then begin
     print,'Vector variables to transform for WREG'
     asknum,'nvector',nvector,doask
     if nvector gt 0 then begin
        vectors=intarr(nvector)
        read,PROMPT='Indices of first components in w? ',vectors
     endif
  endif

end

;===========================================================================
pro do_my_transform,ifile,variables,x,w,xreg,wreg,usereg

  ;; this transformation is useful for plotting RCM files in the
  ;; equatorial plane

  common debug_param & on_error, onerror

  ;; copy x and y coordinates from w to x
  x = w(*,*,0:1)

  ;; change coordinate names
  variables(0:1) = ['x','y']

  ;; use x and w (not xreg and wreg)
  usereg = 0

  ;; Another example: convert 2nd and 3rd indexes from radian to degrees in 3D
  ;; check first if the conversion has been done already
  ;if max(x(*,*,*,1:2)) lt 10.0 then $
  ;  x(*,*,*,1:2) = x(*,*,*,1:2)*180/!pi

end

;===========================================================================
pro do_transform,ifile

  ;; transfrom x, w eithr in-place or to xreg, wreg

  common debug_param & on_error, onerror


  common ask_param, doask
  common file_head              ; gencoord
  common plot_data
  common transform_param

  usereg = (transform eq 'unpolar') or $
           (gencoord and (transform eq 'polar' or transform eq 'regular' $
                          or transform eq 'sphere'))

  if dotransform eq 'n' and (not usereg or keyword_set(wreg)) then return

  if transform eq 'my' or transform eq 'm' then $
     do_my_transform,ifile,variables,x,w,xreg,wreg,usereg $
  else if usereg then case transform of
     'regular': make_regular_grid
     'polar'  :begin
        make_polar_grid
        variables(0:1)=['r','phi']
     end
     'sphere' :begin
        make_sphere_grid
        variables(0:2)=['r','theta','phi']
     end
     'unpolar':begin
        make_unpolar_grid
        variables(0:1)=['x','y']
     end
     else     :print,'Unknown value for transform:',transform
  endcase

end

;===========================================================================
pro read_limits

   common debug_param & on_error, onerror
   
   common ask_param, doask
   common plotfunc_param

   if n_elements(fmax) ne nfunc then fmax=dblarr(nfunc) else fmax=double(fmax)
   if n_elements(fmin) ne nfunc then fmin=dblarr(nfunc) else fmin=double(fmin)

   ; check if there is any function for which autorange is 'y'
   noautorange=1
   for ifunc=0,nfunc-1 do noautorange=noautorange and autoranges(ifunc) eq 'n'

   if(noautorange)then begin
      for ifunc=0,nfunc-1 do begin
         f_min=fmin(ifunc)
         f_max=fmax(ifunc)
         asknum,'Min value for '+funcs(ifunc),f_min,doask
         asknum,'Max value for '+funcs(ifunc),f_max,doask
         fmin(ifunc)=f_min
         fmax(ifunc)=f_max
      endfor
   endif

end
;===========================================================================
pro get_limits,first

  common debug_param & on_error, onerror

  common ask_param, doask
  common plot_data
  common plot_param
  common plotfunc_param
  common transform_param
  common plot_store, f, f1, f2

  for ifunc=0,nfunc-1 do begin
     if autoranges(ifunc) eq 'n' then begin
        if first then begin
           f_min=fmin(ifunc)
           f_max=fmax(ifunc)
           asknum,'Min value for '+funcs(ifunc),f_min,doask
           asknum,'Max value for '+funcs(ifunc),f_max,doask
           fmin(ifunc)=f_min
           fmax(ifunc)=f_max
        endif
     endif else begin
        if usereg then $
           get_func,ifunc,xreg,wreg $
        else $
           get_func,ifunc,x,w

        f_max=max(f)
        f_min=min(f)
        if first then begin
           fmax(ifunc)=f_max
           fmin(ifunc)=f_min
        endif else begin
           if f_max gt fmax(ifunc) then fmax(ifunc)=f_max
           if f_min lt fmin(ifunc) then fmin(ifunc)=f_min
        endelse
     endelse
  endfor

end

;==============================================================================
function get_grid_data,x,y,data,nxreg,xreglimits,triangles,wregpad

  common debug_param & on_error, onerror

  return, griddata(x,y,data,$
                   dimension=nxreg,$
                   start=xreglimits(0:1),$
                   delta=[(xreglimits(2)-xreglimits(0))/(nxreg(0)-1),  $
                          (xreglimits(3)-xreglimits(1))/(nxreg(1)-1)] ,$
                   triangles=triangles,$
;        method='NearestNeighbor',$
;        method='Linear',$
                   method='InverseDistance',$
                   smoothing=0.5,$
                   max_per_sector=4, $
                   missing=wregpad $
                  )
end

;===========================================================================
pro make_regular_grid

;
;    Regularize grid and interpolate "w" via triangulate and trigrid.
;    The original "w" data is interpolated into "wreg", for points outside
;    the convex hull of the irregular grid the "wregpad(nw)" array is used.
;
;    If "x_old" and "x" or "nxreg_old" and "nxreg" are different
;    a triangulization is done first and a regular coordinate array
;    "xreg" is created. The size of the "xreg" array is in "nxreg(2)",
;    "xreglimits(4)" gives the limits. The triangles are saved in "triangles".
;
;    "q" can be interpolated from the irregular grid to the regular one by:
;
;    qreg(*,*)=trigrid(x(*,*,0),x(*,*,1),q,triangles,[0,0],xreglimits)

  common debug_param & on_error, onerror

  common plot_data
  common transform_param
  common file_head              ;  nw

  ;; Floating underflow is not a real error, the message is suppressed
  err=check_math(1,1)

  xx=x(*,*,0)
  yy=x(*,*,1)

  ;; Test distribution. If you uncomment the next lines you can
  ;; take a look at the different "shock wave" representation
  ;; on your grid for the 0-th variable (usually rho)
  ;; for i=0L,n_elements(xx)-1 do $
  ;;   if abs(xx(i))-0.2*abs(yy(i)) gt 0.01 then $
  ;;       w(i,*,0)=2. else w(i,*,0)=1.

  ;; Check if nxreg==nxreg_old and xreglimits==xreglimits_old and x==x_old
  newx=1
  nrectan=0
  if symmtri ne 1 and symmtri ne 2 then $
     if n_elements(nxreg_old) eq n_elements(nxreg) then $
        if max(abs(nxreg_old-nxreg)) eq 0 then $
           if n_elements(xreglimits) eq n_elements(xreglimits_old) then $
              if max(abs(xreglimits-xreglimits_old)) eq 0 then $
                 if n_elements(x_old) eq n_elements(x) then $
                    if max(abs(x_old-x)) eq 0 then newx=0

  if xreglimits(0) eq xreglimits(2) then begin
     xreglimits(0)=min(xx) & xreglimits(2)=max(xx)
  endif
  if xreglimits(1) eq xreglimits(3) then begin
     xreglimits(1)=min(yy) & xreglimits(3)=max(yy)
  endif

  if newx then begin
     print,'Triangulating...'
     x_old=x
     nxreg_old=nxreg
     xreglimits_old=xreglimits

     triangulate,float(xx),float(yy),triangles

                                ; calculate conjugate triangulation and rectangles if required
     if symmtri eq 1 or symmtri eq 2 then $
        symm_triangles,xx,yy,triangles,$
                       triangconj,ntriangles,rectangles,nrectan

     xreg=dblarr(nxreg(0),nxreg(1),2)
     dx=(xreglimits(2)-xreglimits(0))/(nxreg(0)-1)
     dy=(xreglimits(3)-xreglimits(1))/(nxreg(1)-1)
     for i=0,nxreg(1)-1 do xreg(*,i,0)=dx*indgen(nxreg(0))+xreglimits(0)
     for i=0,nxreg(0)-1 do xreg(i,*,1)=dy*indgen(nxreg(1))+xreglimits(1)

  endif
  if n_elements(wregpad) ne nw then begin
     wregpad=dblarr(nw)
     for iw=0,nw-1 do begin
        wmax=max(w(*,*,iw))
        wmin=min(w(*,*,iw))
        if wmax*wmin lt 0 then wregpad(iw)=0 $
        else                   wregpad(iw)=wmin-0.1*(wmax-wmin)
     endfor
  endif

  wreg=dblarr(nxreg(0),nxreg(1),nw)

  case 1 of

     symmtri eq 3: for iw=0,nw-1 do $
        wreg(*,*,iw)=get_grid_data(xx, yy, reform(w(*,*,iw)), nxreg,xreglimits,$
                                   triangles, wregpad(iw))

     symmtri eq 0 or (symmtri lt 3 and nrectan eq 0): for iw=0,nw-1 do $
        wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
                             [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

     symmtri eq 1 and nrectan gt 0: $
        fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
                      triangles,ntriangles,rectangles

     symmtri eq 2 and nrectan gt 0: $
        average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
                          triangles,triangconj

  endcase

  err=check_math(0,0)
                                ;Floating underflow is not a real error, the message is suppressed
  if err ne 32 and err ne 0 then print,'Math error in regulargrid:',err

end

;==============================================================================
pro symm_triangles,xx,yy,triangles,$
                   triangconj,ntriangles,rectangles,nrectan

  common debug_param & on_error, onerror


  ntriangles=n_elements(triangles(0,*))
  print,'Triangulation includes ',ntriangles, '  triangles'
  print,'Checking triangulation ...'

  npoints=n_elements(xx)

  dist=dblarr(npoints-1)
  for i=0L,npoints-2 do $
     dist(i)=(xx(i+1)-xx(i))^2+(yy(i+1)-yy(i))^2
  dist2=min(dist)
  rectangles=lonarr(3,ntriangles)
  ;; Structure of the rectangles array:
  ;; If(rectangles(0,i)=1 then the Ith triangle from the triangles array
  ;; is rectangular one
  tmp1=lonarr(3) & nrec_tri=0
  for i=0L,ntriangles-1 do begin
     if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
            (xx(triangles(1,i))-xx(triangles(2,i)))+$
            (yy(triangles(0,i))-yy(triangles(1,i)))*$
            (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
     then begin
        rectangles(0,i)=1
        tmp1(0)=triangles(1,i)
        if xx(triangles(0,i)) lt xx(triangles(2,i)) then begin
           tmp1(1)=triangles(0,i)
           tmp1(2)=triangles(2,i)
        endif else begin
           tmp1(1)=triangles(2,i)
           tmp1(2)=triangles(0,i)
        endelse
        for j=0,2 do triangles(j,i)=tmp1(j)
     endif

     if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
            (xx(triangles(0,i))-xx(triangles(2,i)))+$
            (yy(triangles(0,i))-yy(triangles(1,i)))*$
            (yy(triangles(0,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
     then begin
        rectangles(0,i)=1
        tmp1(0)=triangles(0,i)
        if xx(triangles(1,i)) lt xx(triangles(2,i)) then begin
           tmp1(1)=triangles(1,i)
           tmp1(2)=triangles(2,i)
        endif else begin
           tmp1(1)=triangles(2,i)
           tmp1(2)=triangles(1,i)
        endelse
        for j=0,2 do triangles(j,i)=tmp1(j)
     endif

     if abs((xx(triangles(0,i))-xx(triangles(2,i)))*$
            (xx(triangles(1,i))-xx(triangles(2,i)))+$
            (yy(triangles(0,i))-yy(triangles(2,i)))*$
            (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
     then begin
        rectangles(0,i)=1
        tmp1(0)=triangles(2,i)
        if xx(triangles(0,i)) lt xx(triangles(1,i)) then begin
           tmp1(1)=triangles(0,i)
           tmp1(2)=triangles(1,i)
        endif else begin
           tmp1(1)=triangles(1,i)
           tmp1(2)=triangles(0,i)
        endelse
        for j=0,2 do triangles(j,i)=tmp1(j)
     endif
  endfor
                                ;Rectangles(1,i) is not equal to zero if the ith rectangular triandgle
                                ;has a common long side with the jth rectangular triangle. In this case
                                ;rectangles(2,i)=j
  nrectan=0L
  for i=0L,ntriangles-1 do begin
     if rectangles(0,i) gt 0 then begin
        nrec_tri=nrec_tri+1
        if rectangles(1,i) eq 0 then begin
           for j=i+1L,ntriangles-1 do begin
              if rectangles(0,j) gt 0 then $
                 if triangles(1,i) eq triangles(1,j) then $
                    if triangles(2,i) eq triangles(2,j) then begin
                 nrectan=nrectan+1
                 rectangles(1,i)=1
                 rectangles(2,i)=j
                 rectangles(1,j)=1
                 rectangles(2,j)=i
                 goto,out
              endif
           endfor
           out:
        endif
     endif
  endfor

  if nrectan ne 0  then begin
     print,'Among    ',nrec_tri, '  rectangular triangles'
     print,'there are',nrectan, '   pairs which have common circumcircle'
     tmp2=lonarr(4)
     ndiag1=0
     ndiag2=0
     triangconj=lonarr(3,ntriangles)
     for i=0L,ntriangles-1 do begin
        if rectangles(1,i) eq 1 then begin
           if rectangles(2,i) gt i then begin
              for j=0,2 do tmp2(j)=triangles(j,i)
              tmp2(3)=triangles(0,rectangles(2,i))
              if yy(tmp2(1)) lt yy(tmp2(2)) then ndiag1=ndiag1+1 else $
                 ndiag2=ndiag2+1
              triangconj(0,i)=tmp2(1)
              triangconj(1,i)=tmp2(0)
              triangconj(2,i)=tmp2(3)
              triangconj(0,rectangles(2,i))=tmp2(2)
              triangconj(1,rectangles(2,i))=tmp2(0)
              triangconj(2,rectangles(2,i))=tmp2(3)
           endif
        endif else for j=0,2 do triangconj(j,i)=triangles(j,i)
     endfor
     print,' Among them ',ndiag1, ' are formed by the triangles,'
     print,' having the common side which is oriented as /////'
     print,' and ',ndiag2, ' have the triangles common side,'
     print,' oriented as \\\\\'
     print,' Calculating the conjugated triangulation ...'
  endif

end

;==============================================================================
pro fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
                  triangles,ntriangles,rectangles

  common debug_param & on_error, onerror

  tmp2=lonarr(4)

  for iw=0,nw-1 do begin
     data=reform(w(*,*,iw))
     print,'Calculating the fitting triangulization for iw=',iw
     for i=0L,ntriangles-1 do begin
        if rectangles(1,i) eq 1 then begin
           if rectangles(2,i) gt i then begin
              tmp2(0:2)=triangles(0:2,i)
              tmp2(3)  =triangles(0,rectangles(2,i))
              if abs(data(tmp2(0))-data(tmp2(3))) lt $
                 abs(data(tmp2(1))-data(tmp2(2))) then begin
                 triangles(0,i)=tmp2(1)
                 triangles(1,i)=tmp2(0)
                 triangles(2,i)=tmp2(3)
                 triangles(0,rectangles(2,i))=tmp2(2)
                 triangles(1,rectangles(2,i))=tmp2(0)
                 triangles(2,rectangles(2,i))=tmp2(3)
              endif
           endif
        endif
     endfor
     wreg(*,*,iw)=trigrid(xx,yy,data,triangles, $
                          [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))
  endfor
  print,'Using fitted triangulation'
end

;==============================================================================
pro average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
                      triangles,triangconj

  common debug_param & on_error, onerror


  wconj=dblarr(nxreg(0),nxreg(1))

  for iw=0,nw-1 do begin
     ;; Calculate wreg with original triangulation
     wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
                          [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1), $
                          missing=wregpad(iw))

     ;; Calculate wconj with conjugated triangulation
     wconj = trigrid(xx,yy,w(*,*,iw),triangconj, $
                     [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1), $
                     missing=wregpad(iw))

     wreg(*,*,iw) = 0.5*(wreg(*,*,iw) + wconj)
  endfor
  print,'Using averaged conjugated triangulation'

end
;===========================================================================
pro make_polar_grid
;
;  Transform 2D or 3D grid and vector variables 
;  from Cartesian to cylindrical components

  common debug_param & on_error, onerror

  common file_head  ; ndim

  case ndim of
      2: make_polar_grid2
      3: make_polar_grid3
      else: begin
          print,'polargid works for 2D and 3D arrays only'
          retall
      end
  endcase
end
;===========================================================================
pro make_polar_grid2
;
;  Transform 2D grid and vector variables 
;  from x and y to radial and phi components

  common debug_param & on_error, onerror

  common plot_data
  common vector_param

  xreg=x
  xreg(*,*,0)=sqrt(x(*,*,0)^2+x(*,*,1)^2)
  xreg(*,*,1)=atan(x(*,*,1),x(*,*,0))
  phi=xreg(*,*,1)
  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)+w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)= -w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor

  ; Remove 2*pi jumps from phi (only works on a regular grid so check ndim)
  sz=size(phi)
  ndim = sz(0)
  if ndim gt 1 then begin
     nx2=sz(2)
     pi2=8*atan(1)
     for ix2=1,nx2-1 do while phi(1,ix2-1) gt phi(1,ix2) do $
        phi(*,ix2)=phi(*,ix2)+pi2
  endif

  xreg(*,*,1)=phi
end

;===========================================================================
pro make_polar_grid3
;
;    Transform 3D vector variables from x,y,z to radial,phi,z components

  common debug_param & on_error, onerror

  common plot_data
  common vector_param

  xreg=x
  xreg(*,*,*,0)=sqrt(x(*,*,*,0)^2+x(*,*,*,1)^2)
  xreg(*,*,*,1)=atan(x(*,*,*,1),x(*,*,*,0))
  phi=xreg(*,*,*,1)
  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,*,ivx)=  w(*,*,*,ivx)*cos(phi)+w(*,*,*,ivy)*sin(phi)
     wreg(*,*,*,ivy)= -w(*,*,*,ivx)*sin(phi)+w(*,*,*,ivy)*cos(phi)
  endfor

  ;Remove 2*pi jumps from phi
  pi2=8*atan(1) & sz=size(phi) & nx2=sz(2)
  for ix2=1,nx2-1 do while phi(1,ix2-1) gt phi(1,ix2) do $
     phi(*,ix2)=phi(*,ix2)+pi2

  xreg(*,*,*,1)=phi
end

;===========================================================================
pro make_sphere_grid
;
;    Transform vector variables from x,y,z to radial,phi,z components

  common debug_param & on_error, onerror

  common plot_data
  common vector_param

  xreg=x
  xreg(*,*,*,0)=sqrt(x(*,*,*,0)^2+x(*,*,*,1)^2+x(*,*,*,2)^2)
  xreg(*,*,*,2)=-atan(x(*,*,*,2),x(*,*,*,0))
  xreg(*,*,*,1)=atan(x(*,*,*,1),sqrt(x(*,*,*,0)^2+x(*,*,*,2)^2))
  phi=xreg(*,*,*,2)
  theta=xreg(*,*,*,1)
  wreg=w
  sinphi=sin(phi)
  cosphi=cos(phi)
  sintheta=sin(theta)
  costheta=cos(theta)
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     ivz=ivy+1
     wreg(*,*,*,ivx)=(w(*,*,*,ivx)*cosphi-w(*,*,*,ivz)*sinphi)*costheta $
                     +w(*,*,*,ivy)*sintheta
     wreg(*,*,*,ivz)=-w(*,*,*,ivx)*sinphi-w(*,*,*,ivz)*cosphi
     wreg(*,*,*,ivy)=(-w(*,*,*,ivx)*cosphi+w(*,*,*,ivz)*sinphi)*sintheta $
                     +w(*,*,*,ivy)*costheta
  endfor

  ;; Remove 2*pi jumps from phi
  pi=4*atan(1) & pi2=2*pi & sz=size(phi) & nx2=sz(2) & nx3=sz(3)
  for ix3=1,nx3-1 do while phi(1,1,ix3-1) gt phi(1,1,ix3) do $
     phi(*,*,ix3)=phi(*,*,ix3)+pi2

  ;Remove turn over from theta
  for ix2=1,nx2-1 do $
  if theta(1,ix2-1,1) ge theta(1,ix2,1) then begin
     if theta(1,ix2,1) lt 0 then $
          theta(*,ix2-1,*)=-pi-theta(*,ix2-1,*) $
     else $
          theta(*,ix2,*)=pi-theta(*,ix2,*)
  endif

  xreg(*,*,*,2)=phi
  xreg(*,*,*,1)=theta
end

;===========================================================================
pro make_unpolar_grid
;
; Transform 2D and 3D vector variables cylindrical to Cartesian components

  common debug_param & on_error, onerror

  common file_head ; ndim

  siz = size(x)
  ndim = siz(0)-1
  case ndim of
      2: make_unpolar_grid2
      3: make_unpolar_grid3
      else: begin
          print,'unpolargid works for 2D and 3D arrays only'
          retall
      end
  endcase

end
;===========================================================================
pro make_unpolar_grid2
;
;    Transform 2D grid (and vector variables)
;    from radial and phi to x and y components

  common debug_param & on_error, onerror

  common plot_data
  common vector_param

  xreg=x
  phi=x(*,*,1)

  phirange = max(phi) - min(phi)

  ; If phi is in local time or degrees, change it to radians
  if      abs(phirange - 24.0) lt 0.1 then phi=phi*!pi/12 $
  else if phirange gt 6.3             then phi=phi*!pi/180

  xreg(*,*,0)=x(*,*,0)*cos(phi)
  xreg(*,*,1)=x(*,*,0)*sin(phi)

  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)-w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)=  w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor
end

;===========================================================================
pro make_unpolar_grid3
;
;  Transform 3D grid (and vector variables)
;  from radial, phi,z to x, y, z components

  common debug_param & on_error, onerror

  common plot_data
  common vector_param

  xreg=x
  phi=x(*,*,*,1)

  ; If phi is in degrees, change it to radians
  if max(abs(phi)) gt 20. then phi=phi*!pi/180

  xreg(*,*,*,0)=x(*,*,*,0)*cos(phi)
  xreg(*,*,*,1)=x(*,*,*,0)*sin(phi)

  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,*,ivx)=  w(*,*,*,ivx)*cos(phi)-w(*,*,*,ivy)*sin(phi)
     wreg(*,*,*,ivy)=  w(*,*,*,ivx)*sin(phi)+w(*,*,*,ivy)*cos(phi)
  endfor
end

;===========================================================================
pro getaxes,ndim,x,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
;===========================================================================
  common debug_param & on_error, onerror

  case ndim of
     1: xx=x
     2: begin
        xx=x(*,*,0)
        yy=x(*,*,1)
     end
     3: begin
        xx=x(*,*,*,0)
        yy=x(*,*,*,1)
        zz=x(*,*,*,2)
     end
  endcase

  if keyword_set(cut0) then begin
     xx=xx(cut0)
     if ndim gt 1 then yy=yy(cut0)
     if ndim gt 2 then zz=zz(cut0)
  endif

  !x.title="!5"+strupcase(variables(0))
  if plotdim gt 1 then !y.title = "!5"+strupcase(variables(1))
  if plotdim gt 2 then !z.title = "!5"+strupcase(variables(2))

  if ndim eq 3 and plotdim eq 2 then begin
     siz=size(cut)
     case 1 of
        siz(0) eq 2: rSlice=zz(0)
        siz(1) eq 1: rSlice=xx(0)
        siz(2) eq 1: rSlice=yy(0)
     endcase
     print,'Normal coordinate of 2D slice:',rSlice
  endif else        rSlice=0.0

; Cut with fixed X value?
  siz=size(cut)
; in 2D
  if siz(0) eq 2 and siz(1) eq 1 then begin
     xx=yy
     !x.title=variables(1)
  endif
; in 3D
  if siz(0) eq 3 then begin
     case 1 of
        plotdim eq 1: begin
           xx=zz
           !x.title=variables(2)
        end
        siz(1) eq 1: begin
           xx=yy
           yy=zz
           !x.title="!5"+variables(1)
           !y.title="!5"+variables(2)
        end
        siz(2) eq 1: begin
           yy=zz
           !y.title="!5"+variables(2)
        end
        else: print,'internal error in getaxes'
     endcase
  endif

end

;===========================================================================
pro set_units, type, distunit=distunit, Mion=Mion, Melectron=Melectron

  ;; If type is given as 
  ;; 'SI', 'CGS', 'NORMALIZED', 'PIC', 'PLANETARY', or 'SOLAR',
  ;; set typeunit = type otherwise try to guess from the fileheader.

  ;; Based on typeunit set units for distance (xSI), time (tSI), 
  ;; density (rhoSI), pressure (pSI), magnetic field (bSI) 
  ;; and current density (jSI) in SI units.
  ;; Distance unit (rplanet or rstar), ion and electron mass in amu 
  ;; can be set with optional distunit, Mion and Melectron.

  ;; Also calculate convenient constants ti0, cs0 ... for typical formulas.
  ;; See file "defaults" for definitions and usage

  common debug_param & on_error, onerror

  common file_head
  common phys_units
  common phys_convert
  common phys_const
  common plot_param ; for rbody

  if keyword_set(type) then $
     typeunit = strupcase(type) $
  else if fixunits then $
     return $
  else if strpos(headline, 'PIC') ge 0 then $
     typeunit = 'PIC' $
  else if strpos(headline,' AU ') ge 0 then $
     typeunit = 'OUTERHELIO' $
  else if strpos(headline, 'kg/m3') ge 0 or strpos(headline,' m/s') ge 0 then $
     typeunit = 'SI' $
  else if strpos(headline,'PLANETARY') ge 0 or strpos(headline,' nPa ') ge 0 or strpos(headline,' nT ') ge 0 then $
     typeunit = 'PLANETARY' $
  else if strpos(headline,' dyne') ge 0 or strpos(headline,' G') ge 0 then $
     typeunit = 'SOLAR' $
  else $
     typeunit = 'NORMALIZED'

  case typeunit of
     'SI': begin
        xSI   = 1.0             ; m
        tSI   = 1.0             ; s
        rhoSI = 1.0             ; kg/m^3
        uSI   = 1.0             ; m/s
        pSI   = 1.0             ; Pa
        bSI   = 1.0             ; T
        jSI   = 1.0             ; A/m^2
     end
     'CGS': begin
        xSI   = 0.01            ; cm
        tSI   = 1.0             ; s
        rhoSI = 1000.0          ; g/cm^3
        uSI   = 0.01            ; cm/s
        pSI   = 0.1             ; dyne/cm^2
        bSI   = 1.0e-4          ; G
        jSI   = 10*cSI          ; Fr/s/cm^2
     end
     'PIC': begin
        ;; Normalized PIC units
        xSI   = 1.0             ; cm
        tSI   = 1.0             ; s
        rhoSI = 1.0             ; g/cm^3
        uSI   = 1.0             ; cm/s
        pSI   = 1.0             ; dyne/cm^2
        bSI   = 1.0             ; G
        jSI   = 1.0             ; Fr/s/cm^2
        c0    = 1.0             ; speed of light always 1 for iPIC3D
     end
     'NORMALIZED': begin
        xSI   = 1.0             ; distance unit in SI
        tSI   = 1.0             ; time unit in SI
        rhoSI = 1.0             ; density unit in SI
        uSI   = 1.0             ; velocity unit in SI
        pSI   = 1.0             ; pressure unit in SI
        bSI   = sqrt(mu0SI)     ; magnetic unit in SI
        jSI   = 1/sqrt(mu0SI)   ; current unit in SI
        c0    = 1.0             ; speed of light (for Boris correction)
     end
     'PLANETARY': begin
        xSI   = ReSi            ; Earth radius (default planet)
        tSI   = 1.0             ; s
        rhoSI = mpSI*1e6        ; mp/cm^3
        uSI   = 1e3             ; km/s
        pSI   = 1e-9            ; nPa
        bSI   = 1e-9            ; nT
        jSI   = 1e-6            ; muA/m^2
        c0    = cSI/uSI         ; speed of light in velocity units
     end
     'OUTERHELIO': begin
        xSI   = AuSI            ; AU
        tSI   = 1.0             ; s
        rhoSI = mpSI*1e6        ; mp/cm^3
        uSI   = 1e3             ; km/s
        pSI   = 1e-1            ; dyne/cm^2
        bSI   = 1e-9            ; nT
        jSI   = 1e-6            ; muA/m^2
        c0    = cSI/uSI         ; speed of light in velocity units
     end
     'SOLAR': begin
        xSI   = RsSI            ; radius of the Sun
        tSI   = 1.0             ; s
        rhoSI = 1e3             ; g/cm^3
        uSI   = 1e3             ; km/s
        pSI   = 1e-1            ; dyne/cm^2
        bSI   = 1e-4            ; G
        jSI   = 1e-6            ; muA/m^2
        c0    = cSI/uSI         ; speed of light in velocity units
     end
     else: begin 
        print, 'ERROR in set_units, invalid typeunit=', typeunit
        retall
     end
  endcase

;; Overwrite values if given by eqpar

  gammae = -1.0

  for ieqpar = 0, neqpar-1 do begin
     value = eqpar(ieqpar)
     case variables(ndim+nw+ieqpar) of
        'xSI'   : xSI    = value
        'tSI'   : tSI    = value
        'uSI'   : uSI    = value
        'rhoSI' : rhoSI  = value
        'mi'    : mi     = value
        'm1'    : mi     = value
        'me'    : me     = value
        'qi'    : qi     = value
        'q1'    : qi     = value
        'qe'    : qe     = value
        'g'     : gamma  = value
        'g1'    : gamma  = value
        'g'     : gamma  = value
        'ge'    : gammae = value
        'c'     : clight = value
        'clight': clight = value
        'r'     : rbody  = value
        'rbody' : rbody  = value
        else:
     endcase
  endfor

  if gammae eq -1 then gammae = gamma

  ;; Overwrite distance unit if given as an argument
  if keyword_set(distunit) then xSI = distunit

  ;; Overwrite ion and electron masses if given as an argument
  if keyword_set(Mion)      then Mi = Mion
  if keyword_set(Melectron) then Me = Melectron

  ;; Calculate convenient conversion factors
  if typeunit eq 'NORMALIZED' then begin
     ti0  = 1.0/Mi              ; T      = p/rho*Mi           = ti0*p/rho
     cs0  = 1.0                 ; cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
     mu0A = 1.0                 ; vA     = sqrt(b/rho)        = sqrt(bb/mu0A/rho)
     mu0  = 1.0                 ; beta   = p/(bb/2)           = p/(bb/(2*mu0))
     uH0  = Mi                  ; uH     = j/rho*Mi           = uH0*j/rho
     op0  = 1.0/Mi              ; omegap = sqrt(rho)/Mi       = op0*sqrt(rho)
     oc0  = 1.0/Mi              ; omegac = b/Mi               = oc0*b
     rg0  = Mi                  ; rg = sqrt(p/rho)/b          = rg0*sqrt(p/rho)/b
     di0  = c0*Mi               ; di = c0/sqrt(rho)*Mi        = di0/sqrt(rho)
     ld0  = Mi                  ; ld = sqrt(p)/(rho*c0)*Mi    = ld0*sqrt(p)/rho
  endif else if typeunit eq 'PIC' then begin
     ti0  = 1.0/Mi              ; T      = p/rho*Mi           = ti0*p/rho        
     cs0  = 1.0                 ; cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)   
     mu0A = 4*!pi               ; vA     = sqrt(b/(4*!pi*rho))= sqrt(bb/mu0A/rho)
     mu0  = 4*!pi               ; beta   = p/(bb/(8*!pi))     = p/(bb/(2*mu0))   
     uH0  = Mi                  ; uH     = j/rho*Mi           = uH0*j/rho        
     op0  = sqrt(4*!pi)/Mi      ; omegap = sqrt(4*!pi*rho)/Mi = op0*sqrt(rho)    
     oc0  = 1.0/Mi              ; omegac = b/Mi               = oc0*b            
     rg0  = Mi                  ; rg = sqrt(p/rho)/b          = rg0*sqrt(p/rho)/b
     di0  = 1.0/sqrt(4*!pi)     ; di = 1/sqrt(4*!pi*rho)*Mi   = di0/sqrt(rho)
     ld0  = 1.0/sqrt(4*!pi)     ; ld = sqrt(p/(4*!pi))/rho*Mi = ld0*sqrt(p)/rho
  endif else begin
     qom  = eSI/(Mi*mpSI) & moq = 1/qom
     ti0  = mpSI/kbSI*pSI/rhoSI*Mi          ; T[K]=p/(nk) = ti0*p/rho
     cs0  = pSI/rhoSI/uSI^2                 ; cs          = sqrt(gs*p/rho)
     mu0A = uSI^2*mu0SI*rhoSI*bSI^(-2)      ; vA          = sqrt(bb/(mu0A*rho))
     mu0  = mu0SI*pSI*bSI^(-2)              ; beta        = p/(bb/(2*mu0))
     uH0  = moq*jSI/rhoSI/uSI               ; uH=j/(ne)   = uH0*j/rho
     op0  = qom*sqrt(rhoSI/e0SI)*tSI        ; omegap      = op0*sqrt(rho)
     oc0  = qom*bSI*tSI                     ; omegac      = oc0*b
     rg0  = moq*sqrt(pSI/rhoSI)/bSI/xSI     ; rg          = rg0*sqrt(p/rho)/b
     di0  = cSI/(op0/tSI)/xSI               ; di=c/omegap = di0/sqrt(rho)
     ld0  = moq*sqrt(pSI)/rhoSI/xSI         ; ld          = ld0*sqrt(p)/rho
  endelse

end

;===========================================================================
pro show_units
;===========================================================================

  common debug_param & on_error, onerror

  common file_head
  common phys_units
  common phys_convert

  print,'headline=', strtrim(headline,2)
  print,'fixunits=', fixunits
  print,'typeunit=', typeunit
  print,'xSI     =', xSI
  print,'tSI     =', tSI
  print,'rhoSI   =', rhoSI
  print,'uSI     =', uSI
  print,'pSI     =', pSI
  print,'bSI     =', bSI
  print,'jSI     =', jSI
  print,'Mi      =', Mi
  print,'Me      =', Me
  print,'Qi      =', Qi
  print,'Qe      =', Qe
  print,'gamma   =', gamma
  print,'gammae  =', gammae
  print,'clight  =', clight
  print,'ti0     =', ti0
  print,'cs0     =', cs0
  print,'mu0A    =', mu0A
  print,'mu0     =', mu0
  print,'c0      =', c0
  print,'uH0     =', uH0
  print,'op0     =', op0
  print,'oc0     =', oc0
  print,'rg0     =', rg0
  print,'di0     =', di0
  print,'ld0     =', ld0

end
;===========================================================================
pro get_func,ifunc,x,w
;===========================================================================
  common debug_param & on_error, onerror

  common plotfunc_param
  common plot_param
  common plot_store

  func1 = funcs1(ifunc)
  func2 = funcs2(ifunc)

  f1 = funcdef(x,w,func1)

  if n_elements(cut0) gt 1 then f1 = f1(cut0)

  if func2 eq '' then f=f1 else begin

     f2 = funcdef(x,w,func2)

     if n_elements(cut0) gt 1 then f2 = f2(cut0)

                                ; Calculate f
     f = sqrt(f1^2 + f2^2)

  endelse

end

;===========================================================================
pro plot_func
;===========================================================================
  common debug_param & on_error, onerror

  common transform_param, usereg
  common plot_data
  common plot_param
  common plotfunc_param
  common plot_store
  common file_head

  ;; Get grid dimensions and set irr=1 
  ;; if it is an irregular grid

  if keyword_set(cut) then siz = size(cut)  $
  else if usereg then      siz = size(xreg) $
  else                     siz = size(x)
  n1=siz(1)
  if plotdim eq 1 then begin
     n2=1
     irr=0
  endif else begin
     n2=siz(2)
     irr=n2 eq 1
  endelse

  if irr and axistype eq 'cells' then begin
     print,'Irregular grid, axistype must be set to coord'
     axistype='coord'
  endif

                                ; Save global values that will be overwritten
  xtitleorig = !x.title
  ytitleorig = !y.title

  if axistype eq 'coord' then begin
     if usereg then $
        getaxes,ndim,xreg,xx,yy,zz,cut,cut0,rSlice,plotdim,variables $
     else $
        getaxes,ndim,x   ,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
  endif

  if xtitleorig ne '' then !x.title=xtitleorig
  if plotdim gt 1 and ytitleorig ne '' then !y.title=ytitleorig

  if !x.range[0] ne !x.range[1] then xrange=!x.range else $
     xrange=[min(xx),max(xx)]
  if plotdim gt 1 then begin
     if !y.range[0] ne !y.range[1] then yrange=!y.range else $
        yrange=[min(yy),max(yy)]
  endif

  ;; Calculate plot spacing from number of subplots per page (ppp) and charsize
  if !p.charsize eq 0.0 then !p.charsize=1.0
  ppp   = multix*multiy
  spacex = float(!d.x_ch_size)/float(!d.x_size)*plot_spacex*!p.charsize
  spacey = float(!d.y_ch_size)/float(!d.y_size)*plot_spacey*!p.charsize
  set_space, ppp, spacex, spacey, sizes, nx = multix, ny = multiy

  ;; Store x and y titles and tick names
  newytitle  = 0 ;;; ytitleorig ne !y.title
  xtitle     = !x.title
  ytitle     = !y.title
  xtickname  = !x.tickname
  ytickname  = !y.tickname

  for ifunc=0,nfunc-1 do begin

     plotmod=plotmodes(ifunc)

                                ; stream2 --> stream
     i=strpos(plotmod,'stream2')
     if i ge 0 then plotmod=strmid(plotmod,0,i+6)+strmid(plotmod,i+7)

                                ; contour --> cont
     i=strpos(plotmod,'contour')
     if i ge 0 then plotmod=strmid(plotmod,0,i+4)+strmid(plotmod,i+7)

     i=strpos(plotmod,'white')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+5)
        white=1
     endif else white=0

     i=strpos(plotmod,'noaxis')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+6)
        noaxis=4
     endif else noaxis = 0

     i=strpos(plotmod,'grid')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
        showgrid=1
     endif else showgrid=0

     i=strpos(plotmod,'irr')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
        irr=1
     endif

     i=strpos(plotmod,'mesh')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
        showgrid=1
        if irr then showmesh=0 else showmesh=1
     endif else showmesh=0

     i=strpos(plotmod,'body')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
        showbody=1
     endif else showbody=0

     i=strpos(plotmod,'map')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
        showmap=1
     endif else showmap=0

     i=strpos(plotmod,'usa')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
        showusa=1
     endif else showusa=0

     i=strpos(plotmod,'fill')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
        fill=1
     endif else fill=0

     i=strpos(plotmod,'bar')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
        showbar=1
        if strpos(plotmod,'stream') ge 0 then showbar=0
        fill=1
     endif else showbar=0

     i=strpos(plotmod,'label')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+5)
        label=1
     endif else label=0

     i=strpos(plotmod,'log')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
        logarithm=1
     endif else logarithm=0

     i=strpos(plotmod,'max')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
        usemax=1
     endif else usemax=0

     i=strpos(plotmod,'mean')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
        usemean=1
     endif else usemean=0

                                ; check if this plot requires a special color table #ctNNN
     i = strpos(plotmod, '#ct')
     if i ge 0 then begin
        color = strmid(plotmod,i+3)
        plotmod = strmid(plotmod,0,i)
        loadct_extra, color
     endif

                                ; check if this plot requires a special color #cNNN
     i = strpos(plotmod, '#c')
     if i ge 0 then begin
        color = strmid(plotmod,i+2)
        plotmod = strmid(plotmod,0,i)+strmid(plotmod,i+5)
        !p.color = color
     endif

     !p.title=plottitles(ifunc)
     if !p.title eq 'default' then !p.title=funcs(ifunc)
     !p.title="!5"+!p.title

     if ifunc lt nfunc-1 then begin
        if strpos(plotmodes(ifunc+1),'over') gt 0 $
           and not keyword_set(timetitle) then begin
           nexttitle = plottitles(ifunc+1)
           if nexttitle eq 'default' then nexttitle=funcs(ifunc+1)
           if nexttitle ne ' ' then begin
              if strpos(plotmodes(ifunc+1),'log') gt 0 then $
                 nexttitle = 'log '+nexttitle
              !p.title += ' and '+nexttitle
           endif
        endif
     endif

     i=strpos(plotmod,'over')
     if i ge 0 then begin
        plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
        !p.multi(0) = !p.multi(0)+1
        if !p.multi(0) ge !p.multi(1)*!p.multi(2) then !p.multi(0)=0
        !p.title = ''
     endif

     ;; Check if the angular unit of phi is given
     angleunit = -1.0
     i=strpos(plotmod,'polardeg')
     if i ge 0 then begin
        angleunit = !dtor
        plotmod=strmid(plotmod,0,i+5)+strmid(plotmod,i+8)
     endif
     i=strpos(plotmod,'polarhour')
     if i ge 0 then begin
        angleunit = !pi/12
        plotmod=strmid(plotmod,0,i+5)+strmid(plotmod,i+9)
     endif
     i=strpos(plotmod,'polarrad')
     if i ge 0 then begin
        angleunit = 1.0
        plotmod=strmid(plotmod,0,i+5)+strmid(plotmod,i+8)
     endif

     ;; Calculate the next p.multi(0) explicitly
     if !p.multi(0) gt 0 then multi0=!p.multi(0)-1 $
     else multi0=!p.multi(1)*!p.multi(2)-1

     ;; Calculate subplot position indices
     if !p.multi(4) then begin
        ;; columnwise
        plotix=multix-1-multi0/multiy
        plotiy=multi0 mod multiy
     endif else begin
        ;; rowwise
        plotix=multix-1-(multi0 mod multix)
        plotiy=multi0/multix
     endelse

     if plotix eq 0 then newytitle = 0

     if plotmod ne 'shade' and plotmod ne 'surface' then begin

        ;; obtain position for flat plotmodes
        set_position, sizes, plotix, multiy-1-plotiy, pos, /rect

        ;; shrink in X direction if there is a colorbar among the plotmodes
        if strpos(plotmode,'bar') gt 0 then $
           pos(2) = pos(2) - (pos(2) - pos(0))*0.15

        ;; shrink in X direction for the Y axis of plot
        if strmid(plotmod,0,4) eq 'plot' and multix gt 1 then $
           pos(0) = pos(0) + (pos(2) - pos(0))*0.15

        if keyword_set(fixaspect) and strmid(plotmod,0,4) ne 'plot' then begin

           if plotmod eq 'polar' then $
              aspectx=1 $
           else if abs(fixaspect) ne 1 then $
              aspectx = abs(fixaspect) $
           else begin
              if axistype eq 'cells' then begin
                 width  = n1-1.0
                 height = n2-1.0
              endif else begin
                 width  = abs(xrange(1) - xrange(0))
                 height = abs(yrange(1) - yrange(0))
              endelse

              aspectx = width/height
           endelse

           aspectpos = (pos(2)-pos(0))/(pos(3)-pos(1)) $
                       *float(!d.x_size)/float(!d.y_size)

           aspectratio = aspectpos/aspectx

           if aspectratio gt 1 then begin
              if fixaspect gt 0 then begin
                 posmid=(pos(2)+pos(0))/2.
                 posdif=(pos(2)-pos(0))/2.
                 pos(0)=posmid - posdif/aspectratio
                 pos(2)=posmid + posdif/aspectratio
              endif else $
                 pos(0:2:2) = pos(0:2:2)/aspectratio + 0.5*(1-1/aspectratio)
           endif else begin
              if fixaspect gt 0 then begin
                 posmid=(pos(3)+pos(1))/2.
                 posdif=(pos(3)-pos(1))/2.
                 pos(1)=posmid - posdif*aspectratio
                 pos(3)=posmid + posdif*aspectratio
              endif else $
                 pos(1:3:2) = pos(1:3:2)*aspectratio + 0.5*(1-aspectratio)
           endelse

        endif

        ;; Omit X axis if unneeded
        if plotiy gt 0 and not showxtitle then begin
           if not showxaxis then !x.tickname = strarr(60)+' '
           !x.title = ' '
        endif

        ;; Omit Y axis if unneeded
        if plotix gt 0 and not showytitle $
           and strmid(plotmod,0,4) ne 'plot' then begin
           if not showyaxis then !y.tickname = strarr(60)+' '
           !y.title = ' '
        endif

        !p.position = pos

     endif     

     if usereg then get_func,ifunc,xreg,wreg $
     else           get_func,ifunc,x,w

                                ; calculate running max or mean
     if nplotstore gt 0 then begin
                                ; allocate plotstore and timesture arrays
        if n_elements(plotstore) ne $
           n_elements(f)*nplotstore*nfunc*nfilestore then $
              plotstore = fltarr(n_elements(f),nplotstore,nfunc,nfilestore)

        if n_elements(timestore) ne nplotstore*nfilestore then $
           timestore = fltarr(nplotstore,nfilestore)

                                ; store plot function (list of points) and simulation time
        plotstore(*,iplotstore,ifunc,ifilestore) = f
        timestore(iplotstore,ifilestore) = time

                                ; jump to next storage index (cyclic)
        if ifunc eq nfunc-1 then iplotstore = (iplotstore + 1) mod nplotstore

                                ; calculate maximum over stored functions
        if usemax then $
           for i = 0, nplotstore-1 do f = f > plotstore(*,i,ifunc,ifilestore)

                                ; calculate mean of stored functions
        if usemean then begin
           f = 0*f
           for i = 0, nplotstore-1 do f = f + plotstore(*,i,ifunc,ifilestore)
           f = f/nplotstore
        endif
     endif

     f_min = fmin(ifunc)
     f_max = fmax(ifunc)

     if logarithm and f_min gt 0 and min(f) gt 0 then begin
        f     = alog10(f)
        f_min = alog10(f_min)
        f_max = alog10(f_max)
        if plottitles(ifunc) eq 'default' and !p.title ne '' $
        then !p.title = 'log '+!p.title
     endif

     if usemax and plottitles(ifunc) eq 'default' then $
        !p.title = 'max '+!p.title

     if usemean and plottitles(ifunc) eq 'default' then $
        !p.title = 'mean '+!p.title

     if f_max eq f_min then begin
        f_max = f_max + 1
        f_min = f_min - 1
     endif

     if strmid(plotmod,0,4) eq 'plot' then $
        if nfunc gt ppp                then lstyle=ifunc/ppp $
        else if keyword_set(linestyle) then lstyle=linestyle $
        else                                lstyle=!p.linestyle

     ;; Skip minimum ad maximum levels
     if plotmod eq 'cont' or plotmod eq 'polar' then begin
        if fill then nlevel=colorlevel else nlevel=contourlevel
        levels=(findgen(nlevel+2)-1)/(nlevel-1)*(f_max-f_min)+f_min
     endif

     ;; figure out the units of angle in 2nd coordinate if not already set
     if plotmod eq 'polar' and angleunit lt 0 then begin
        if max(yy)-min(yy) gt 300 then $
           angleunit = !dtor $  ; degrees
        else if max(yy)-min(yy) gt 20 then $
           angleunit = !pi/12 $ ; local time
        else $
           angleunit = 1.0      ; radians
     endif

     if plotmod eq 'tv' then begin
        ;; Calculate plotting position and size
        tvplotx=pos(0)*!d.x_size
        tvploty=pos(1)*!d.y_size
        tvsizex=(pos(2)-pos(0))*!d.x_size
        tvsizey=(pos(3)-pos(1))*!d.y_size
                                ; recalculate f for tv mode
        if !d.name eq 'PS' then tvf=congrid(f,200,200) $
        else                    tvf=congrid(f,tvsizex,tvsizey)

        tvf=bytscl(tvf,MIN=f_min,MAX=f_max,TOP=!D.TABLE_SIZE-3)+1
     endif

     if showbar then plot_color_bar, $
        [pos(2)+0.005, pos(1), pos(2)+0.025, pos(3)], [f_min,f_max]

     case axistype of
        'cells': case plotmod of
           'cont': contour,f>f_min,LEVELS=levels,$
                           FILL=fill,FOLLOW=label,$
                           XSTYLE=noaxis+1,YSTYLE=noaxis+1,/NOERASE
           'plot':plot,f,YRANGE=[f_min,f_max],$
                       XSTYLE=noaxis+18,ystyle=18,LINE=lstyle,/NOERASE
           'plot_io':plot_io,f,YRANGE=[f_min,f_max],$
                             XSTYLE=noaxis+18,ystyle=18,LINE=lstyle,/NOERASE
           'shade'    :begin
              shade_surf,f>f_min,ZRANGE=[f_min,f_max],$
                         XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                         ZSTYLE=noaxis+18,AX=viewanglex,AZ=viewanglez,/NOERASE
              if showgrid then $
                 surface,f>f_min,ZRANGE=[f_min,f_max],$
                         XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                         ZSTYLE=noaxis+18,AX=viewanglex,AZ=viewanglez,/NOERASE
           end
           'surface'  :surface,f>f_min,ZRANGE=[f_min,f_max],$
                               XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                               ZSTYLE=noaxis+18,AX=viewanglex,AZ=viewanglez,/NOERASE
           'tv'       :begin
              tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
              contour,f,XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                      /NODATA,/NOERASE
           end
           'vel'      :vector, f1, f2, NVECS=velvector, MAXVAL=f_max, $
                              DYNAMIC=velspeed, SEED=velseed, X0=velpos, $
                              /NOERASE, WHITE=white, $
                              XSTYLE=noaxis+1, YSTYLE=noaxis+1
           'vector'   :vector, f1, f2, NVECS=velvector, MAXVAL=f_max, $
                              DYNAMIC=velspeed, SEED=velseed, X0=velpos, $
                              /NOERASE,WHITE=white, $
                              XSTYLE=noaxis+1, YSTYLE=noaxis+1
           'stream'   :streamline, f1, f2, NVECS=velvector, SEED=velseed, $
                                  X0=velpos, /NOERASE, WHITE=white, $
                                  XSTYLE=noaxis+1, YSTYLE=noaxis+1
           'velovect' :velovect,f1,f2,/NOERASE
           'ovelovect':velovect,f1,f2,/NOERASE,$
                                XRANGE=[0,n_elements(f1(*,0))-1], $
                                YRANGE=[0,n_elements(f1(0,*))-1]
        endcase
        'coord': case plotmod of
           'cont'     :if irr then begin
              if not keyword_set(tri) then $
                 triangulate,float(xx),float(yy),tri
              contour,f>f_min,xx,yy,$
                      FOLLOW=label, FILL=fill, TRIANGULATION=tri, $
                      LEVELS=levels,XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                      /NOERASE
           endif else $
              contour,f>f_min,xx,yy,$
                      FOLLOW=label, FILL=fill, LEVELS=levels, $
                      XSTYLE=noaxis+1,YSTYLE=noaxis+1,/NOERASE
           'polar'    :polar_contour,f>f_min,yy*angleunit,xx,$
                                     FOLLOW=label, FILL=fill, LEVELS=levels,$
                                     XSTYLE=noaxis+1,YSTYLE=noaxis+1,/dither, $
                                     /NOERASE
           'plot'     :plot,xx,f,YRANGE=[f_min,f_max],$
                            XSTYLE=noaxis+1,YSTYLE=noaxis+3,$
                            LINE=lstyle,/NOERASE
           'plot_io'  :plot_io,xx,f,YRANGE=[f_min,f_max],$
                               XSTYLE=noaxis+1,YSTYLE=noaxis+3,$
                               LINE=lstyle,/NOERASE
           'plot_oi'  :plot_oi,xx,f,YRANGE=[f_min,f_max],$
                               XSTYLE=noaxis+1,YSTYLE=noaxis+3,$
                               LINE=lstyle,/NOERASE
           'plot_oo'  :plot_oo,xx,f,YRANGE=[f_min,f_max],$
                               XSTYLE=noaxis+1,YSTYLE=noaxis+3,$
                               LINE=lstyle,/NOERASE
           'shade'    :if irr then begin
              shade_surf_irr,f>f_min,xx,yy,AX=viewanglex,AZ=viewanglez
              shade_surf,f>f_min,xx,yy,AX=viewanglex,AZ=viewanglez,/NODATA,/NOERASE
           endif else begin
              shade_surf,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                         XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                         ZSTYLE=noaxis+18,AX=viewanglex,AZ=viewanglez,/NOERASE
              if showgrid then $
                 surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                         XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                         ZSTYLE=noaxis+18,$
                         AX=viewanglex,AZ=viewanglez,/NOERASE
           endelse
           'surface'  :surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                               XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                               ZSTYLE=noaxis+18,$
                               AX=viewanglex,AZ=viewanglez,/NOERASE
           'tv'       :begin
              tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
              contour,f,xx,yy,$
                      XSTYLE=noaxis+1,YSTYLE=noaxis+1,$
                      /NODATA,/NOERASE
           end
           'vel'      :vector, f1, f2, xx, yy, XXOLD=velx, YYOLD=vely, $
                               TRIANGLES=veltri, NVECS=velvector,MAXVAL=f_max,$
                               DYNAMIC=velspeed, SEED=velseed, X0=velpos,$
                               /NOERASE, WHITE=white, $
                               XSTYLE=noaxis+1, YSTYLE=noaxis+1
           'vector'   :vector, f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                               TRIANGLES=veltri,NVECS=velvector,MAXVAL=f_max,$
                               DYNAMIC=velspeed, SEED=velseed, X0=velpos,$
                               /NOERASE, WHITE=white, $
                               XSTYLE=noaxis+1, YSTYLE=noaxis+1
           'arrow'    :vector, f1/f, f2/f, xx, yy, XXOLD=velx, YYOLD=vely, $
                               TRIANGLES=veltri, NVECS=velvector, MAXVAL=1.5,$
                               DYNAMIC=velspeed, SEED=velseed, X0=velpos,$
                               /NOERASE, WHITE=white, $
                               XSTYLE=noaxis+1, YSTYLE=noaxis+1
           'stream'   :streamline, f1, f2, xx, yy, XXOLD=velx, YYOLD=vely,$
                                   TRIANGLES=veltri, NVECS=velvector,$
                                   SEED=velseed, X0=velpos,$
                                   /NOERASE, WHITE=white, $
                                   XSTYLE=noaxis+1, YSTYLE=noaxis+1
           'velovect' :velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE
           'ovelovect':velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE,$
                                XRANGE=xrange, YRANGE=yrange
        endcase
        else:print,'Unknown axistype:',axistype
     endcase

     if showbody and axistype eq 'coord' then $
        if rBody gt abs(rSlice) then begin
        rBody = float(rBody)    ; make sure it's not an integer
        theta = findgen(37)*!pi*2.0/36.0
        rBodySlice=sqrt(rBody^2-rSlice^2)

        polyfill, rBodySlice*cos(theta), rBodySlice*sin(theta),color = 0, $
                  noclip=0
                                ; redraw box in case the body is at the edge
        if(plotmod ne 'polar')then $
           plot,xx,yy,XSTYLE=noaxis+1,YSTYLE=noaxis+1,/NODATA,/NOERASE
     endif

     if showmap or showusa then $
        map_set, 0.0, 180., $
                 /cylindrical, /continent, usa=showusa, /noborder, /noerase, $
                 limit=[yrange(0),xrange(0),yrange(1),xrange(1)]

     !p.title = ''

     if showgrid and plotdim eq 2 and plotmod ne 'surface'    $
        and plotmod ne 'shade' then begin

        if white then tvlct,bytarr(256,3)+255 ; change to all white color table

        if(plotmod eq 'polar')then                                       $
           plot_grid,xx,yy*angleunit,lines=showmesh,xstyle=5,ystyle=5,/polar $
        else if keyword_set(cut) then                                    $
           plot_grid,xx,yy,lines=showmesh,xstyle=5,ystyle=5               $
        else begin
           plot_grid,x,lines=showmesh,xstyle=5,ystyle=5,$
                     xrange=xrange,yrange=yrange
        endelse

                                ; restore colors
        if white then tvlct, r_curr, g_curr, b_curr

     endif

     !p.multi(0) = multi0
     !p.position = 0
     !x.title    = xtitle
     !x.tickname = xtickname
     !y.title    = ytitle
     !y.tickname = ytickname
  endfor

  !p.position = 0
  !x.title    = xtitleorig
  !y.title    = ytitleorig

end
;===========================================================================
pro putbottom,multix,multiy,ix,iy,info,nx,it,time

  common debug_param & on_error, onerror

  t    = round(time)
  sec  = t mod 60
  t    = t/60
  min  = t mod 60
  t    = t/60
  hour = t mod 24
  t    = t/24
  day  = t mod 365
  year = t/365

  if year gt 0 then $
     stime = string(year,day,hour,format='(i3,"y",i3.3,"d",i2.2,"h")') $
  else if day gt 0 then $
     stime = string(day,hour,min,format='(i4,"d",i2.2,"h",i2.2,"m")') $
  else if hour gt 0 then $
     stime = string(hour,min,sec,format='(i4,"h",i2.2,"m",i2.2,"s")') $
  else $
     stime = string(time, format='(g12.5)')

  if n_elements(nx) eq 1 then snx = string(nx, format='(i6)')
  if n_elements(nx) eq 2 then snx = string(nx, format='(i6,",",i4)')
  if n_elements(nx) eq 3 then snx = string(nx, format='(i6,",",i4,",",i4)')

  sit = string(it, format='(i8)')

  siz = size(info)

  if siz(1) eq 7 then begin
     if not execute("s = "+info) then begin
        print,'Error in putbottom: cannot evaluate info=',info
        return
     endif
  endif else begin
     if info lt 1 then return
     s = ''
     if info gt 2 then s = 'nx=' + snx + ', '
     if info gt 1 then s = s + 'it=' + sit + ', '
     s = s + 'time=' + stime
  endelse

  xyouts, 5+(ix*!d.x_size)/multix, 8+(iy*!d.y_size)/multiy, /DEV, s

;;xyouts,5,3000,/DEV,info

;info='time ='+string(time/(17.24*3600),format='(f6.2)')+' Uday'
;xyouts,3000+(0.9*ix*!d.x_size)/multix,3000+(0.82*(iy+1)*!d.y_size)/multiy,/DEV,info,charsize=1.2

;info='time ='+string(time/60,format='(i4)')+' min'
;xyouts,!p.position[0]-0.1*(!p.position(2)-!p.position(0)),$
;       !p.position[1]-0.1*(!p.position(3)-!p.position(1)),$
;       info,/NORM,charsize=0.8

end
;===========================================================================
pro putheader,multix,multiy,ix,iy,ninfo,headline,nx

  common debug_param & on_error, onerror

  if ninfo lt 1 then return
  info=strtrim(headline,2)
  if ninfo gt 1 then info=info+' (nx='+string(nx,format='(i6,2(i4))')+')'
  xyouts,5+(ix*!d.x_size)/multix,-12+((iy+1)*!d.y_size)/multiy,/DEV,info

end
;===========================================================================
function diff1,a,x
;
; Take derivative of "a" with respect to "x" (if present)
; using 2nd order centered differencing
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  ndim=siz(0)
  if ndim ne 1 then begin
     print,'Function diff1 is intended for 1D arrays only'
     retall
  endif

  n  = n_elements(a)
  nx = n_elements(x)

  if nx ne 0 and nx ne n then begin
     print,'Error in diff1, arrays sizes differ: nx, n=', nx, n
     retall
  endif

  dadx = a

  if nx eq 0 then dadx(1:n-2) = (a(2:n-1) - a(0:n-3))/2 $
  else            dadx(1:n-2) = (a(2:n-1) - a(0:n-3))/(x(2:n-1) - x(0:n-3))

; fill in boundaries
  dadx(0)   = dadx(1)
  dadx(n-1) = dadx(n-2)

  return,dadx
end
;===========================================================================
function laplace,a,x,y,z
;
; Take Laplace of "a" in 1 or 2D with respect to "x" (and "y") if present.
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  ndim=siz(0)
  if ndim eq 3 then return,laplace3(a,x,y,z)
  if ndim eq 2 then return,laplace2(a,x,y)
  if ndim ne 1 then begin
     print,'Function laplace is intended for 1, 2 and 3D arrays only'
     retall
  endif

  n  = n_elements(a)
  nx = n_elements(x)

  if nx ne 0 and nx ne n then begin
     print,'Error in laplace, array sizes differ: nx, n=', nx, n
     retall
  endif

  d2adx2 = a

  if nx eq 0 then d2adx2(1:n-2) = a(2:n-1) - 2*a(1:n-2) + a(0:n-3) $
  else            d2adx2(1:n-2) = $
     ( (a(2:n-1) - a(1:n-2))/(x(2:n-1) - x(1:n-2)) $
       - (a(1:n-2) - a(0:n-3))/(x(1:n-2) - x(0:n-3)) ) $
     / (0.5*(x(2:n-1) - x(0:n-3)))

; fill in boundaries
  d2adx2(0)   = d2adx2(1)
  d2adx2(n-1) = d2adx2(n-2)

  return,d2adx2
end
;===========================================================================
function laplace2,a,x,y
;
; Take Laplace of "a" in 2D with respect to "x" and "y" (if present).
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  ndim=siz(0)

  n1 = siz(1)
  n2 = siz(2)
  n  = n1*n2
  nx = n_elements(x)
  ny = n_elements(y)

  if (nx ne 0 or ny ne 0) and (nx ne n or ny ne n) then begin
     print,'Error in laplace2, array sizes differ: n, nx, ny=', n, nx, ny
     retall
  endif

  d2adx2 = a

  if nx eq 0 then d2adx2(1:n1-2,1:n2-2) = $
     a(2:n1-1,1:n2-2) - 2*a(1:n1-2,1:n2-2) + a(0:n1-3,1:n2-2) + $
     a(1:n1-2,2:n2-1) - 2*a(1:n1-2,1:n2-2) + a(1:n1-2,0:n2-3)   $
  else            d2adx2(1:n1-2,1:n2-2) = $
     ( (a(2:n1-1,1:n2-2) - a(1:n1-2,1:n2-2))/ $
       (x(2:n1-1,1:n2-2) - x(1:n1-2,1:n2-2))- $
       (a(1:n1-2,1:n2-2) - a(0:n1-3,1:n2-2))/ $
       (x(1:n1-2,1:n2-2) - x(0:n1-3,1:n2-2)) $
     ) / (0.5*(x(2:n1-1,1:n2-2) - x(0:n1-3,1:n2-2))) + $
     ( (a(1:n1-2,2:n2-1) - a(1:n1-2,1:n2-1))/ $
       (y(1:n1-2,2:n2-1) - y(1:n1-2,1:n2-1))- $
       (a(1:n1-2,1:n2-2) - a(1:n1-2,0:n2-3))/ $
       (y(1:n1-2,1:n2-2) - y(1:n1-2,0:n2-3)) $
     ) / (0.5*(y(1:n1-2,2:n2-1) - y(1:n1-2,0:n2-3)))

; fill in boundaries
  d2adx2(0   ,1:n2-1) = d2adx2(1   ,1:n2-1)
  d2adx2(n1-1,1:n2-1) = d2adx2(n1-2,1:n2-1)
  d2adx2(*,0)         = d2adx2(*,1)
  d2adx2(*,n2-1)      = d2adx2(*,n2-2)

  return,d2adx2
end
;===========================================================================
function laplace3,a,x,y,z
;
; Take Laplace of "a" in 2D with respect to "x", "y" and "z" (if present).
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  ndim=siz(0)

  n1 = siz(1)
  n2 = siz(2)
  n3 = siz(3)
  n  = n1*n2*n3
  nx = n_elements(x)
  ny = n_elements(y)
  nz = n_elements(z)

  if (nx+ny+nz ne 0) and (nx ne n or ny ne n or nz ne n) then begin
     print,'Error in laplace3, array sizes differ: n, nx, ny, nz=', $
           n, nx, ny, nz
     retall
  endif

  d2adx2 = a

  if nx eq 0 then d2adx2(1:n1-2,1:n2-2,1:n3-2) = $
     a(2:n1-1,1:n2-2,1:n3-2) + a(0:n1-3,1:n2-2,1:n3-2) + $
     a(1:n1-2,2:n2-1,1:n3-2) + a(1:n1-2,0:n2-3,1:n3-2) + $
     a(1:n1-2,1:n2-2,2:n3-1) + a(1:n1-2,1:n2-2,0:n3-3) - $
     6*a(1:n1-2,1:n2-2,1:n3-2)                           $
  else            d2adx2(1:n1-2,1:n2-2,1:n3-2) = $
     ( (a(2:n1-1,1:n2-2,1:n3-2) - a(1:n1-2,1:n2-2,1:n3-2))/ $
       (x(2:n1-1,1:n2-2,1:n3-2) - x(1:n1-2,1:n2-2,1:n3-2))- $
       (a(1:n1-2,1:n2-2,1:n3-2) - a(0:n1-3,1:n2-2,1:n3-2))/ $
       (x(1:n1-2,1:n2-2,1:n3-2) - x(0:n1-3,1:n2-2,1:n3-2)) $
     ) / (0.5*(x(2:n1-1,1:n2-2,1:n3-2) - x(0:n1-3,1:n2-2,1:n3-2))) + $
     ( (a(1:n1-2,2:n2-1,1:n3-2) - a(1:n1-2,1:n2-1,1:n3-2))/ $
       (y(1:n1-2,2:n2-1,1:n3-2) - y(1:n1-2,1:n2-1,1:n3-2))- $
       (a(1:n1-2,1:n2-2,1:n3-2) - a(1:n1-2,0:n2-3,1:n3-2))/ $
       (y(1:n1-2,1:n2-2,1:n3-2) - y(1:n1-2,0:n2-3,1:n3-2)) $
     ) / (0.5*(y(1:n1-2,2:n2-1,1:n3-2) - y(1:n1-2,0:n2-3,1:n3-2))) + $
     ( (a(1:n1-2,1:n2-2,2:n3-1) - a(1:n1-2,1:n2-2,1:n3-2))/ $
       (z(1:n1-2,1:n2-2,2:n3-1) - z(1:n1-2,1:n2-2,1:n3-2))- $
       (a(1:n1-2,1:n2-2,1:n3-2) - a(1:n1-2,1:n2-2,0:n3-3))/ $
       (z(1:n1-2,1:n2-2,1:n3-2) - z(1:n1-2,1:n2-2,0:n3-3)) $
     ) / (0.5*(z(1:n1-2,1:n2-2,2:n3-1) - z(1:n1-2,1:n2-2,0:n3-3)))

; fill in boundaries
  d2adx2(0   ,1:n2-1,1:n3-2) = d2adx2(1   ,1:n2-1,1:n3-2)
  d2adx2(n1-1,1:n2-1,1:n3-2) = d2adx2(n1-2,1:n2-1,1:n3-2)
  d2adx2(*   ,0     ,1:n3-2) = d2adx2(*   ,1     ,1:n3-2)
  d2adx2(*   ,n2-1  ,1:n3-2) = d2adx2(*   ,n2-2  ,1:n3-2)
  d2adx2(*   ,*     ,0     ) = d2adx2(*   ,*     ,1     )
  d2adx2(*   ,*     ,n3-1  ) = d2adx2(*   ,*     ,n3-2  )

  return,d2adx2
end
;===========================================================================
function diff2,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 2nd order centered differencing
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  ndim=siz(0)
  if ndim ne 2 and ndim ne 3 then begin
     print,'Function diff2 is intended for 2D and 3D arrays only'
     retall
  endif

  if direction lt 1 or direction gt ndim then begin
     print,'Direction=',direction,' should be between 1 and ndim=',ndim,'!'
     retall
  endif

  n1=siz(1)
  n2=siz(2)
  if ndim eq 3 then n3=siz(3)

  if direction eq 1 then begin
     ind1=indgen(n1)
     jnd1=ind1+1
     jnd1(n1-1)=n1
     hnd1=ind1-1
     hnd1(0)=0
     if ndim eq 2 then $
        dadx=(a(jnd1,*)-a(hnd1,*))/(x(jnd1,*)-x(hnd1,*)) $
     else $
        dadx=(a(jnd1,*,*)-a(hnd1,*,*))/(x(jnd1,*,*)-x(hnd1,*,*))
  endif
  if direction eq 2 then begin
     ind2=indgen(n2)
     jnd2=ind2+1
     jnd2(n2-1)=n2
     hnd2=ind2-1
     hnd2(0)=0
     if ndim eq 2 then $
        dadx=(a(*,jnd2)-a(*,hnd2))/(x(*,jnd2)-x(*,hnd2)) $
     else $
        dadx=(a(*,jnd2,*)-a(*,hnd2,*))/(x(*,jnd2,*)-x(*,hnd2,*))
  endif
  if direction eq 3 then begin
     ind3=indgen(n3)
     jnd3=ind3+1
     jnd3(n3-1)=n3
     hnd3=ind3-1
     hnd3(0)=0
     dadx=(a(*,*,jnd3)-a(*,*,hnd3))/(x(*,*,jnd3)-x(*,*,hnd3))
  endif

  return,dadx

end

;===========================================================================

function diff2avg,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 2nd order centered differencing and average in the orthogonal direction.
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  ndim=siz(0)
  if ndim ne 2 then begin
     print,'Function diff2avg is intended for 2D arrays only'
     retall
  endif

  if direction lt 1 or direction gt ndim then begin
     print,'Direction=',direction,' should be between 1 and ndim=',ndim,'!'
     retall
  endif

  n1=siz(1)
  n2=siz(2)

  avg = a
  
  ind1=indgen(n1)
  jnd1=ind1+1
  jnd1(n1-1)=n1-1
  hnd1=ind1-1
  hnd1(0)=0

  ind2=indgen(n2)
  jnd2=ind2+1
  jnd2(n2-1)=n2-1
  hnd2=ind2-1
  hnd2(0)=0

  if direction eq 1 then begin    
     avg = 0.25*a(*,hnd2) + 0.5*a + 0.25*a(*,jnd2)
     dadx = (avg[jnd1,*] - avg[hnd1,*])/(x[jnd1,*] - x[hnd1,*]) 
  endif

  if direction eq 2 then begin
     avg = 0.25*a[hnd1,*] + 0.5*a + 0.25*a[jnd1,*]
     dadx = (avg[*,jnd2] - avg[*,hnd2])/(x[*,jnd2] - x[*,hnd2]) 
  endif
  
  return,dadx

end

;===========================================================================

function diff4,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 4th order centered differencing
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  if siz(0) ne 2 then begin
     print,'Function diff4 is intended for 2D arrays only'
     retall
  endif

  n1=siz(1)
  n2=siz(2)

  dadx=a

  if direction eq 1 then begin
     if n1 lt 5 then begin
        print,'Cannot take 4th order X gradient of grid with less than 5 columns'
        retall
     endif
     dadx(2:n1-3,*)=(a(4:n1-1,*)-8*a(3:n1-2,*)+8*a(1:n1-4,*)-a(0:n1-5,*)) $
                    /(x(3:n1-2,*)-x(1:n1-4,*))/6 ;
     dadx(0,*)   =dadx(2,*)
     dadx(1,*)   =dadx(2,*)
     dadx(n1-2,*)=dadx(n1-3,*)
     dadx(n1-1,*)=dadx(n1-3,*)
  endif
  if direction eq 2 then begin
     if n2 lt 5 then begin
        print,'Cannot take 4th order Y gradient of grid with less than 5 rows'
        retall
     endif
     dadx(*,2:n2-3)=(a(*,4:n2-1)-8*a(*,3:n2-2)+8*a(*,1:n2-4)-a(*,0:n2-5)) $
                    /(x(*,3:n2-2)-x(*,1:n2-4))/6 ;
     dadx(*,0)   =dadx(*,2)
     dadx(*,1)   =dadx(*,2)
     dadx(*,n2-2)=dadx(*,n2-3)
     dadx(*,n2-1)=dadx(*,n2-3)
  endif

  return,dadx

end

;===========================================================================
function diff3,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using IDL's 1D deriv() function
;
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(a)
  if siz(0) ne 2 then begin
     print,'Function diff3 is intended for 2D arrays only'
     retall
  endif

  dadx=a

  if direction eq 1 then for i2=0,siz(2)-1 do dadx(*,i2)=deriv(x(*,i2),a(*,i2))
  if direction eq 2 then for i1=0,siz(1)-1 do dadx(i1,*)=deriv(x(i1,*),a(i1,*))

  return,dadx

end

;===========================================================================
function minmod,a,b
;
; Calculate minmod limited slope of a and b slopes

  common debug_param & on_error, onerror

  ;; get sign of a
  if a gt 0 then s=1 else s=-1

  ;; calculate limited slope
  c = s*max([0,min([abs(a),s*b])])

  return,c
end
;==========================================================================
function symmdiff,direction,a,x,y,report=report,anti=anti
;
;find the symmetric for irregular grid
;
;=========================================================================
  
  common debug_param & on_error, onerror

  ;; If x is not present, the grid is taken to be Cartesian
  if n_elements(x) eq 0 then return, symmdiffreg(direction,a,anti=anti)

  if not keyword_set(report) then report = 0

  diff=a*0

  e = exp(1.0d0)
  sorting = x + e*y
  sortsym = x - e*y
  n = n_elements(x)

  if direction eq 1 then begin
     signx = -1.0 & signy = +1.0
  endif else begin
     signx = +1.0 & signy = -1.0
  endelse

  for i=0,n-1 do begin
     if diff(i) eq 0 then begin
        xi = x(i)
        yi = y(i)   
        sortsymi = sortsym(i)
        j  = n/2
        dj = j
                                ;print,'i,x,y,sortsymi=',i,xi,yi,sortsymi
                                ;print,'n,j,dj=',n,j,dj
        last = 0
        fail = 0
        while abs(x(j)-signx*xi) gt 1e-6 or abs(y(j)-signy*yi) gt 1e-6 do begin
           if last then begin
              fail = 1
              break
           endif
           
                                ;print,'j,dj,sorting=',j,dj,sorting(j)

           if dj eq 1 then last = 1 else dj = (dj+1)/2

           if sortsymi le sorting(j) then j = j - dj else j = j + dj

           if j lt 0 then j = 0
           if j ge n then j = n - 1

        endwhile 

                                ;print,'solution j,x,y,sort=',j,x(j),y(j),sorting(j)

        if fail then begin
           print,'error in symdiff: incorrect number of mirror points'
           print,'i,j,dj,last=',i,j,dj,last
           print,'xi,yi=',xi,yi
           if j ge 0 and j lt n then print,'xj,yj=',x(j),y(j)
           diff(i)=1e30
           retall
        endif else begin
           if keyword_set(anti) then diff(i)= a(i)+a(j) else diff(i)= a(i)-a(j)
           diff(j)=-diff(i)
           if report and abs(diff(i)) gt report then print,"i,x,y,diff=",i,xi,yi,diff(i)
        endelse
     endif
  endfor
  return,diff

end
;===========================================================================
function symmdiffreg,direction,a,anti=anti
;
; Take symmetric difference of "a" with respect to a mirror plane in direction
; "direction"
;
;===========================================================================
  common debug_param & on_error, onerror

  if not keyword_set(report) then report = 0

  siz=size(a)
  dim=siz(0)
  nx=siz(1)

  diff=a

  coef = 1.0
  if keyword_set(anti) then coef = -1.0

  case dim of
     1: for i=0,nx-1 do diff(i) = a(i) - coef*a(nx-1-i)
     2: begin
        ny=siz(2)
        case direction of
           1: for i=0,nx-1 do diff(i,*) = a(i,*)-coef*a(nx-1-i,*)
           2: for i=0,ny-1 do diff(*,i) = a(*,i)-coef*a(*,ny-1-i)
        endcase
     end
     3: begin
        ny=siz(2)
        nz=siz(3)
        case direction of
           1: for i=0,nx-1 do diff(i,*,*) = a(i,*,*) - coef*a(nx-1-i,*,*)
           2: for i=0,ny-1 do diff(*,i,*) = a(*,i,*) - coef*a(*,ny-1-i,*)
           3: for i=0,nz-1 do diff(*,*,i) = a(*,*,i) - coef*a(*,*,nz-1-i)
        endcase
     end
     4: begin
        ny=siz(2)
        nz=siz(3)
        case direction of
           1: for i=0,nx-1 do diff(i,*,*,*)=a(i,*,*,*)-coef*a(nx-1-i,*,*,*)
           2: for i=0,ny-1 do diff(*,i,*,*)=a(*,i,*,*)-coef*a(*,ny-1-i,*,*)
           3: for i=0,nz-1 do diff(*,*,i,*)=a(*,*,i,*)-coef*a(*,*,nz-1-i,*)
        endcase
     end
  endcase

  return,diff
end
;===========================================================================
function filledge,a

; On the edges use copy of closest cells

  common debug_param & on_error, onerror

  siz=size(a)
  n1=siz(1)
  n2=siz(2)

  result=a
  result(0,*)   =result(1,*)
  result(*,0)   =result(*,1)
  result(n1-1,*)=result(n1-2,*)
  result(*,n2-1)=result(*,n2-2)

  return,result
end

;===========================================================================
pro gengrid,name,x,y,xc,yc,vol2,u,v
;
; From cell center coordinates x,y calculate cell corner coordinates xc,yc,
; cell volumes. Check for array sizes of the optional u,v arguments.
; The name of the calling function is shown for error messages.
;===========================================================================
  common debug_param & on_error, onerror

  siz=size(x)
  if siz(0) ne 2 then begin
     print,'Function ',name,' is for 2D arrays only'
     retall
  endif

  n1=siz(1)
  n2=siz(2)

  error=''
  siz=size(y)
  if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd coord'
  if keyword_set(u) then begin
     siz=size(u)
     if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='1st func'
  endif
  if keyword_set(v) then begin
     siz=size(v)
     if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd func'
  endif
  if error ne '' then begin
     print,'In function ',name,' the first argument does not match the ',error
     retall
  endif

; Coordinates for cell corners
  xc=(x(0:n1-2,0:n2-2)+x(0:n1-2,1:n2-1)+x(1:n1-1,0:n2-2)+x(1:n1-1,1:n2-1))/4
  yc=(y(0:n1-2,0:n2-2)+y(0:n1-2,1:n2-1)+y(1:n1-1,0:n2-2)+y(1:n1-1,1:n2-1))/4

; Calculate 2*volume=(diagonal_1 X diagonal_2)
  vol2=dblarr(n1,n2)+1
  vol2(1:n1-2,1:n2-2)= $
     ((xc(1:n1-2,1:n2-2)-xc(0:n1-3,0:n2-3))*(yc(0:n1-3,1:n2-2)-yc(1:n1-2,0:n2-3)) $
      -(yc(1:n1-2,1:n2-2)-yc(0:n1-3,0:n2-3))*(xc(0:n1-3,1:n2-2)-xc(1:n1-2,0:n2-3)))

end

;===========================================================================
function intedge,f,xc
;
; Integrate the neighbouring values of "f" for the four edges described by "xc"
; The size of "f", "xc", and the result are n1*n2, (n1-1)*(n2-1), and n1*n2
; respectively, but only the inner (n1-2)*(n2-2) points are calculated, the
; edge values are 0-s.
;===========================================================================

  common debug_param & on_error, onerror

  siz=size(f)
  n1=siz(1)
  n2=siz(2)

  intf=dblarr(n1,n2)
  intf(1:n1-2,1:n2-2)=-(xc(1:n1-2,1:n2-2)-xc(0:n1-3,1:n2-2))*f(1:n1-2,2:n2-1) $
                      -(xc(1:n1-2,0:n2-3)-xc(1:n1-2,1:n2-2))*f(2:n1-1,1:n2-2) $
                      -(xc(0:n1-3,0:n2-3)-xc(1:n1-2,0:n2-3))*f(1:n1-2,0:n2-3) $
                      -(xc(0:n1-3,1:n2-2)-xc(0:n1-3,0:n2-3))*f(0:n1-3,1:n2-2)

  return,intf

end

;===========================================================================
function intedge_rz,f,rc,zc
;
; Integrate r_edge*f_neighbour*dz for the four cell edges.
; assuming axial symmetry in the ignored direction.
; Only the inner (n1-2)*(n2-2) points are calculated, the edge values are 0-s.
;
;===========================================================================

  common debug_param & on_error, onerror

  siz=size(f)
  n1=siz(1)
  n2=siz(2)

  intf=dblarr(n1,n2)
  intf(1:n1-2,1:n2-2)= $
     -f(1:n1-2,2:n2-1)*(rc(1:n1-2,1:n2-2)+rc(0:n1-3,1:n2-2)) $
     *(zc(1:n1-2,1:n2-2)-zc(0:n1-3,1:n2-2)) $
     -f(2:n1-1,1:n2-2)*(rc(1:n1-2,0:n2-3)+rc(1:n1-2,1:n2-2)) $
     *(zc(1:n1-2,0:n2-3)-zc(1:n1-2,1:n2-2)) $
     -f(1:n1-2,0:n2-3)*(rc(0:n1-3,0:n2-3)+rc(1:n1-2,0:n2-3)) $
     *(zc(0:n1-3,0:n2-3)-zc(1:n1-2,0:n2-3)) $
     -f(0:n1-3,1:n2-2)*(rc(0:n1-3,1:n2-2)+rc(0:n1-3,0:n2-3)) $
     *(zc(0:n1-3,1:n2-2)-zc(0:n1-3,0:n2-3))

  return,intf

end

;===========================================================================
function grad_2d,idir,f,x,y
;
; Take gradient of "f" in direction "idir" on the "x,y" structured 2D grid.
; Gradient is the contour integral of edge_normal_idir*f_edge_averaged
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Gradient can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

  common debug_param & on_error, onerror

  if n_elements(idir) eq 0 or n_elements(f) eq 0 $
     or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
     print,'Missing arguments in function grad'
     retall
  endif

  gengrid,'grad',x,y,xc,yc,vol2,f

  if idir eq 1 then return,filledge( intedge(f,yc)/vol2) $
  else              return,filledge(-intedge(f,xc)/vol2)

end

;===========================================================================
function grad_rz,idir,f,r,z
;
; Take gradient of "f" in direction "idir" on the "r,z" structured 2D grid
; assuming axial symmetry in the ignored direction.
; Gradient is the contour integral of edge_normal_idir*f*R_edge_averaged
; divided by R*cell_volume - f/r for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels for idir=2, or equals +f/2R for idir=1.
; Gradient can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

  common debug_param & on_error, onerror

  if n_elements(ndir) eq 0 or n_elements(f) eq 0 $
     or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
     print,'Missing arguments in function grad_rz'
     retall
  endif

  gengrid,'grad_rz',r,z,rc,zc,vol2,f

  if idir eq 1 then return,filledge( (intedge_rz(f,rc,zc)/vol2 - f)/2/r ) $
  else              return,filledge( -intedge(f,rc^2)/vol2/2/r)

end

;===========================================================================
function div,u,v,x,y
;
; Take divergence of "u,v" vector with respect to "x,y" on a structured 2D grid
; Divergence is the contour integral of edge_normal.(u,v)_edge_averaged
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Divergence can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

  common debug_param & on_error, onerror

  if n_elements(u) eq 0 or n_elements(v) eq 0 $
     or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
     print,'Missing arguments in function div'
     retall
  endif

  gengrid,'div',x,y,xc,yc,vol2,u,v

  return,filledge((intedge(u,yc)-intedge(v,xc))/vol2)

end

;===========================================================================
function div_rz,u,v,r,z
;
; Take divergence of "u,v" vector with respect to "r,z" on a structured 2D grid
; assuming axial symmetry in the ignored direction.
; Divergence is the contour integral of edge_normal.(u,v)*R_edge_averaged
; divided by R*cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center is simply u/(2R).
; Divergence can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

  common debug_param & on_error, onerror

  if n_elements(u) eq 0 or n_elements(v) eq 0 $
     or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
     print,'Missing arguments in function div_rz'
     retall
  endif

  gengrid,'div_rz',r,z,rc,zc,vol2,u,v

  return,filledge(((intedge_rz(u,rc,zc)-intedge(v,rc^2))/vol2 + u)/2/r)

end

;===========================================================================
function curl,u,v,x,y
;
; Take curl of "u,v" vector with respect to "x,y" on a structured 2D grid.
; Curl is the contour integral of edge_vector.(u,v)_edge_averaged
; divided by cell_volume for each cells. See also comments for div function.
;
;===========================================================================

  common debug_param & on_error, onerror

  if n_elements(u) eq 0 or n_elements(v) eq 0 $
     or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
     print,'Missing arguments in function curl'
     retall
  endif

  gengrid,'curl',x,y,xc,yc,vol2,u,v

  return,filledge((intedge(u,xc)+intedge(v,yc))/vol2)

end

;===========================================================================
function curl_rz,u,v,r,z
;
; Take curl of "u,v" vector with respect to "r,z" on a structured 2D grid
; with axial symmetry in the ignored direction.
; Curl is the contour integral of edge_vector.(u,v)*R_edge_averaged
; divided by R*cell_volume for each cells - v/R.
; See also comments for the div_rz function on edge average and edge cells.
;
;===========================================================================
  common debug_param & on_error, onerror

  if n_elements(u) eq 0 or n_elements(v) eq 0 $
     or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
     print,'Missing arguments in function curl_rz'
     retall
  endif

  gengrid,'curl',r,z,rc,zc,vol2,u,v

  return,filledge(-((intedge_rz(v,rc,zc)+intedge(u,rc^2))/vol2 - v)/2/r)

end

;===========================================================================
function quadruplet,nx,x0,x1,dx,ny,y0,y1,dy,nz,z0,z1,dz,nw,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=quadruplet(100,0,30,2,100,30,40,1)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(quadruplet(100,0,99,4,100,30,69,2,2,0,1,1))
;===========================================================================

  common debug_param & on_error, onerror

  if keyword_set(dx) then begin
     checkdim,1,nx,x0,x1,dx
     all=lindgen(x1+1)
     sub=all(x0:x1)
     ind=sub(where(sub mod dx eq x0 mod dx))
  end
  if keyword_set(dy) then begin
     checkdim,2,ny,y0,y1,dy
     ixs=ind
     all=lindgen(y1+1)
     sub=all(y0:y1)
     iys=sub(where(sub mod dy eq y0 mod dy))
     ind=(ixs # (0*iys+1)) + ((0*ixs+nx) # iys)
  end
  if keyword_set(dz) then begin
     checkdim,3,nz,z0,z1,dz
     ixys=ind
     nxy=long(nx)*long(ny)
     all=lindgen(z1+1)
     sub=all(z0:z1)
     izs=sub(where(sub mod dz eq z0 mod dz))
     ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs))
     for iz=0,n_elements(izs)-1 do ind(*,*,iz)=ixys + izs(iz)*nxy
  end
  if keyword_set(dw) then begin
     checkdim,4,nw,w0,w1,dw
     ixyzs=ind
     nxyz=long(nx)*long(ny)*long(nz)
     all=lindgen(w1+1)
     sub=all(w0:w1)
     iws=sub(where(sub mod dw eq w0 mod dw))
     ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs),n_elements(iws))
     for iw=0,n_elements(iws)-1 do ind(*,*,*,iw)=ixyzs + iws(iw)*nxyz
  end

  return,ind
end

;===========================================================================
function triplet,x0,x1,dx,y0,y1,dy,z0,z1,dz,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=triplet(0,99,2,0,99,2)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(triplet(0,99,4,0,99,4,0,1,1))
;
; Note: the resulting indices are valid for an array of size
;
;      (x1+1)*(y1+1)*(z1+1)*(w1+1)
;
;===========================================================================

  common debug_param & on_error, onerror

  if keyword_set(dw) then $
     return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz,w1+1,w0,w1,dw)

  if keyword_set(dz) then $
     return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz)

  if keyword_set(dy) then $
     return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy)

  if keyword_set(dx) then $
     return,quadruplet(x1+1,x0,x1,dx)

  print,'Error in TRIPLET: All strides are 0!'
  retall

end

;==========================================
function coarsen,a,boxsize,fd=fd
;
; Produce a coarser array from "a" by averaging out cells in a box.
; The box size can be defined by a scalar (n long interval, n*n square,
; or n*n*n cube or n*n*n*n hyper cube) or as an
; array of the same dimension as "a" (n1*n2 rectangle or n1*n2*n3
; brick or n1*n2*n3*n4 hyper brick)
;
; If /fd is set then use finite difference coarsening:
; extract every n1-th element in dimension 1, every n2-th in dim 2 etc.

  common debug_param & on_error, onerror

  if(n_elements(a) eq 0 or n_elements(boxsize) eq 0)then begin
     print,'Calling sequence is: array_co=coarse(array, boxsize, /fd)'
     retall
  endif

  siz=size(a)
  ndim=siz(0)

  if(ndim eq 0 or ndim gt 4)then begin
     print,'coarse requires a 1, 2, 3 or 4D array for the 1st argument'
     retall
  endif
  nx=siz(1:ndim)

  siz=size(box)
  if(siz(0) eq 0)then begin
     n = intarr(ndim) + boxsize
  endif else if siz(0) eq ndim then begin
     n = boxsize
  endif else begin
     print,'boxsize should either be a scalar, or an array '
     print,'of the same dimension as the number of dimensions of the array'
     retall
  endelse

  if keyword_set(fd) then case ndim of
     1: result = a(triplet(0,nx(0)-1,n(0)))
     2: result = a(triplet(0,nx(0)-1,n(0), 0,nx(1)-1,n(1)))
     3: result = a(triplet(0,nx(0)-1,n(0), 0,nx(1)-1,n(1), 0,nx(2)-1,n(2)))
     4: result = a(triplet(0,nx(0)-1,n(0), 0,nx(1)-1,n(1), 0,nx(2)-1,n(2), $
                           0,nx(3)-1,n(3)))
  endcase else case ndim of
     1: begin
        result = dblarr(nx(0)/n(0))
        for ix=0,nx(0)/n(0)-1 do $
           for i=0,n(0)-1 do $
              result(ix)=result(ix) + a(ix*n(0)+i)
        result=result/n(0)
     end
     2: begin
        result = dblarr(nx(0)/n(0),nx(1)/n(1))
        for ix=0,nx(0)/n(0)-1 do $
           for iy=0,nx(1)/n(1)-1 do $
              for i=0,n(0)-1 do $
                 for j=0,n(1)-1 do $
                    result(ix,iy) = result(ix,iy) + a(ix*n(0)+i,iy*n(1)+j)
        result=result/n(0)/n(1)
     end
     3: begin
        result=dblarr(nx(0)/n(0),nx(1)/n(1),nx(2)/n(2))
        for ix=0,nx(0)/n(0)-1 do $
           for iy=0,nx(1)/n(1)-1 do $
              for iz=0,nx(2)/n(2)-1 do $
                 for i=0,n(0)-1 do $
                    for j=0,n(1)-1 do $
                       for k=0,n(2)-1 do $
                          result(ix,iy,iz) = result(ix,iy,iz) $
           + a(ix*n(0)+i,iy*n(1)+j,iz*n(2)+k)
        result = result/n(0)/n(1)/n(2)
     end
     4: begin
        result = dblarr(nx(0)/n(0),nx(1)/n(1),nx(2)/n(2),nx(3)/n(3))
        for ix=0,nx(0)/n(0)-1 do $
           for iy=0,nx(1)/n(1)-1 do $
              for iz=0,nx(2)/n(2)-1 do $
                 for iw=0,nx(3)/n(3)-1 do $
                    for i=0,n(0)-1 do $
                       for j=0,n(1)-1 do $
                          for k=0,n(2)-1 do $
                             for l=0,n(3)-1 do $
                                result(ix,iy,iz,iw) = result(ix,iy,iz,iw) $
           + a(ix*n(0)+i,iy*n(1)+j,iz*n(2)+k,iw*n(3)+l)
        result = result/n(0)/n(1)/n(2)/n(3)
     end
  endcase

  return,result
end

;===========================================================================
pro checkdim,idim,nx,x0,x1,dx

; Check quadruplet for conditions nx>x1>=x0>=0 and dx>0
;===========================================================================
  common debug_param & on_error, onerror

  if nx le 0 then begin
     print,'Size must be positive for dimension',idim
     retall
  endif
  if x1 ge nx then begin
     print,'Maximum index must be less than size for dimension',idim
     retall
  endif
  if x0 lt 0 then begin
     print,'Minimum index must be greater than 0 for dimension',idim
     retall
  endif
  if x0 gt x1 then begin
     print,'Minimum index must be less than maximum index for dimension',idim
     retall
  endif
  if dx le 0 then begin
     print,'Stride must be a positive integer for dimension',idim
     retall
  endif

  return
end

;===================================================================
pro plot_grid,x,y,lines=lines,xstyle=xstyle,ystyle=ystyle,polar=polar,$
              xrange=xrange,yrange=yrange,noorigin=noorigin
;===================================================================

  common debug_param & on_error, onerror

  if not keyword_set(x) then begin
     print,'Usage: plot_grid, x [,y] [,/lines] [,/polar]',$
           ' [,xstyle=3] [,ystyle=1]',$
           '                   [,xrange=[-10,10]], [yrange=[-10,10]]'
     retall
  endif

  xx=reform2(x)
  sizx=size(xx)

  if (n_elements(polar) eq 0) then polar = 0

  if not keyword_set(y) then begin
     case sizx(0) of
        3:begin
           if sizx(3) ne 2 then goto, ERROR1
           yy=xx(*,*,1)
           xx=xx(*,*,0)
        end
        2:begin
           if sizx(2) ne 2 then goto, ERROR1
           yy=xx(*,1)
           xx=xx(*,0)
           lines=0
        end
        else: goto, ERROR1
     endcase
  endif else begin
     yy=reform2(y)
     sizy=size(yy)
     if sizx(0) ne sizy(0)            then goto, ERROR2
     if max(abs(sizx-sizy)) ne 0      then goto, ERROR2
     if sizx(0) ne 2 and sizx(0) ne 1 then goto, ERROR2
     if sizx(0) eq 1 then lines=0
  endelse

  if not keyword_set(xrange) then xrange = [0,0]
  if not keyword_set(yrange) then yrange = [0,0]

  if keyword_set(lines) then begin

     plot, xx, yy, XSTYLE=xstyle, YSTYLE=ystyle, POLAR=polar, $
           XRANGE=xrange, YRANGE=yrange, /NOERASE, /NODATA

     if(keyword_set(noorigin))then begin
        for ix=0,sizx(1)-1 do $
           for iy=0,sizx(2)-2 do $
              if((xx(ix,iy)   ne 0 or yy(ix,iy)   ne 0) and $
                 (xx(ix,iy+1) ne 0 or yy(ix,iy+1) ne 0)) then $
                    oplot,[xx(ix,iy),xx(ix,iy+1)],[yy(ix,iy),yy(ix,iy+1)],POLAR=polar,$
                          psym=0

        for iy=0,sizx(2)-1 do $
           for ix=0,sizx(1)-2 do $
              if((xx(ix,iy)   ne 0 or yy(ix,iy)   ne 0) and $
                 (xx(ix+1,iy) ne 0 or yy(ix+1,iy) ne 0)) then $
                    oplot,[xx(ix,iy),xx(ix+1,iy)],[yy(ix,iy),yy(ix+1,iy)],POLAR=polar,$
                          psym=0

     endif else begin

        for ix=0,sizx(1)-1 do $
           oplot,xx(ix,*),yy(ix,*),POLAR=polar,psym=0
        for iy=0,sizx(2)-1 do $
           oplot,xx(*,iy),yy(*,iy),POLAR=polar,psym=0

     endelse

  endif else begin

     if polar then $
        plot, xx, yy, PSYM=3, SYMSIZE=!p.symsize, $
              XRANGE=xrange, YRANGE=yrange, XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE,$
              /POLAR $
     else $
        plot, xx, yy, PSYM=1, SYMSIZE=!p.symsize, $
              XRANGE=xrange, YRANGE=yrange, XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE
  endelse

  return

ERROR1:
  print,'size(x)=',sizx
  print,'Error: plot_grid,x  requires x(nx,ny,2) array'
  retall

ERROR2:
  print,'size(x)=',sizx,' size(y)=',sizy
  print,'Error: plot_grid,x,y requires x(nx,ny) y(nx,ny) arrays'
  retall


end

;==========================================
pro compare,w0,w1,wnames

; Compare all variables in w0 and w1 by calculating
; relative difference in the 1st norm.
;==========================================

  common debug_param & on_error, onerror

  sizew0=size(w0)
  sizew1=size(w1)

  if sizew0(0) ne sizew1(0) then begin
     print,'w0 and w1 have different dimensions:',sizew0(0),' and ',sizew1(0)
     retall
  endif

  ndim=sizew0(0)-1

  if ndim eq 0 then begin
     ndim=1
     nw=1
  endif else $
     nw=sizew0(ndim+1)

  if max(abs(sizew0(1:ndim)-sizew1(1:ndim))) gt 0 then begin
     print,'w0 and w1 have different sizes:',sizew0(1:ndim),' /= ',sizew1(1:ndim)
     retall
  endif

  if keyword_set(wnames) then $
     print,$
     'var 2*max(|A-B|)/max(|A|+|B|) 2*sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)/2' $
  else $
     print, $
     'ind 2*max(|A-B|)/max(|A|+|B|) 2*sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)/2'

  for iw=0,nw-1 do begin
     case ndim of
        1: begin
           wsum=max(abs(w0(*,iw))+abs(w1(*,iw)))/2
           wdif=max(abs(w0(*,iw)-w1(*,iw)))
           wsum1=total(abs(w0(*,iw))+abs(w1(*,iw)))/2
           wdif1=total(abs(w0(*,iw)-w1(*,iw)))
        end
        2: begin
           wsum=max(abs(w0(*,*,iw))+abs(w1(*,*,iw)))/2
           wdif=max(abs(w0(*,*,iw)-w1(*,*,iw)))
           wsum1=total(abs(w0(*,*,iw))+abs(w1(*,*,iw)))/2
           wdif1=total(abs(w0(*,*,iw)-w1(*,*,iw)))
        end
        3: begin
           wsum=max(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))/2
           wdif=max(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
           wsum1=total(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))/2
           wdif1=total(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
        end
     endcase

     if keyword_set(wnames) then begin
        if wsum eq 0. then print,wnames(iw),' wsum=0' $
        else               print,wnames(iw),wdif/wsum,wdif1/wsum1,wsum
     endif else begin
        if wsum eq 0. then print,iw,' wsum=0' $
        else               print,iw,wdif/wsum,wdif1/wsum1,wsum
     endelse

  endfor
end

;=============================================================================
function rel_error, w1, w2, iws, fd=fd

; Calculate relative errors of w1 with respect to w2.
; The 1 to ndim (=1,2, or 3) indexes are the spatial indexes. 
; The last index is assumed to be the variable index.
; If w1 has more elements than w2, it is regarded as the "reference"
; solution, otherwise w2 is regarded as reference solution.
; The difference between w1 and w2 is normalized by the reference solution.
; The function returns the relative errors averaged over all variables.
; If iws is present then the error is calculated for the variable
; indexes listed in iws. 
;
; If the /fd keyword is set, use finite difference coarsening.

  common debug_param & on_error, onerror

  if n_elements(w1) le n_elements(w2) then begin
     w     = w1
     wref  = w2
  endif else begin
     w     = w2
     wref  = w1
  endelse

  s     = size(w)
  sref  = size(wref)
  ndim  = s(0) - 1
  nw    = s(ndim+1)

  if nw ne sref(ndim+1) then begin
     print,'ERROR in rel_error: w1 and w2 have different number of variables'
     help,w1,w2
     retall
  endif

  if n_elements(iws) gt 0 then begin
     ; keep variable indexes iws only in w, wref
     nw = n_elements(iws)
     case ndim of
        1: begin
           w    = w(*,iws)
           wref = wref(*,iws)
        end
        2: begin
           w    = w(*,*,iws)
           wref = wref(*,*,iws)
        end
        3: begin
           w    = w(*,*,*,iws)
           wref = wref(*,*,*,iws)
        end
     endcase
  endif

  ; Coarsen wref if necessary (only in the spatial indexes)
  nx    = s(1:ndim)
  nxref = sref(1:ndim)
  if max(nxref gt nx) then wref = coarsen(wref, [nxref/nx, 1], fd=fd)

  error = double(0.0)
  for iw = 0, nw-1 do begin
     case ndim of
        1: error = error + (total(abs(w(*,iw)     - wref(*,iw)))) $
                   /        total(abs(wref(*,iw)))
        2: error = error + (total(abs(w(*,*,iw)   - wref(*,*,iw)))) $
                   /        total(abs(wref(*,*,iw)))
        3: error = error + (total(abs(w(*,*,*,iw) - wref(*,*,*,iw)))) $
                   /        total(abs(wref(*,*,*,iw)))
     endcase
  endfor

  return, error/nw

end

;=============================================================================
function rel_errors, w0, w1, w2, w3, w4, w5, ivar=ivar, ratio=ratio, fd=fd

  ; Calculate relative errors for up to 6 arrays w0 ...
  ; The arrays can have different sizes. The finer arrays are coarsened. 
  ; The last array should be the reference solution.
  ;
  ; The ivar keyword lists the variables to be compared. Default is all.
  ; If /ratio keyword is set and there are at least 3 arguments, 
  ;    then print out the ratio of errors. Useful for convergence rate.
  ; If /fd keyword is set, use finite difference style coarsening.
  ;    The default is finite volume style coarsening.


  common debug_param & on_error, onerror

  n = lonarr(6)
  n(0) = n_elements(w0)
  n(1) = n_elements(w1)
  n(2) = n_elements(w2)
  n(3) = n_elements(w3)
  n(4) = n_elements(w4)
  n(5) = n_elements(w5)

  narray = fix(total(n gt 0))
  if narray lt 2 then begin
     print,'ERROR in rel_errors: number of elements of arrays=',n
     print,'There should be at least 2 non-empty arrays!'
     retall
  endif

  nmax = max(n)
  if narray gt 2 and nmax ne n(narray-1) then begin
     print,'ERROR in rel_errors: number of elements=',n
     print,'The last array should have the most elements!'
     retall
  endif

  if narray eq 2 then wref = w1
  if narray eq 3 then wref = w2
  if narray eq 4 then wref = w3
  if narray eq 5 then wref = w4
  if narray eq 6 then wref = w5

  errors = dblarr(narray-1)

  errors(0)                     = rel_error(w0, wref, ivar, fd=fd)
  if narray gt 2 then errors(1) = rel_error(w1, wref, ivar, fd=fd)
  if narray gt 3 then errors(2) = rel_error(w2, wref, ivar, fd=fd)
  if narray gt 4 then errors(3) = rel_error(w3, wref, ivar, fd=fd)
  if narray gt 5 then errors(4) = rel_error(w4, wref, ivar, fd=fd)

  if keyword_set(ratio) and narray gt 2 then $
     print,'ratio=',errors(0:narray-3)/(errors(1:narray-2) > 1d-25)

  return,errors
end

;=============================================================================
pro get_log, source, wlog, wlognames, logtime, timeunit, verbose=verbose

; Read the log data from source. If source is an integer, it is 
; interpreted as a unit number. If it is a string, it is taken as the
; filename. Read the content of the file into wlog and the variable 
; names into wlognames. 
; The optional logtime argument is set to the time in hours.
; If verbose is present set show verbose information.
; If versbose is a string, attach it to 'wlog' in the verbose info.

  common debug_param & on_error, onerror

  if not keyword_set(source) then begin
     print, $
        'Usage: get_log, source, wlog, wlognames [,logtime, timeunit] [,verbose=verbose]'
     help,source,wlog,wlognames
     retall
  endif

  itype = size(source,/type)
  if itype eq 2 or itype eq 3 then begin
     filesource=0
     unit = source
     file = 'unit '+strtrim(string(unit),2)
     stat = fstat(unit)
     if not stat.open then begin
        print,'get_log error: unit is not open'
        retall
     endif
  endif else if itype eq 7 then begin
     filesource=1
     file = source
     unit = 0
     found = 0
     while not found do begin
        unit = unit + 1
        stat = fstat(unit)
        if not stat.open then found = 1
     endwhile
     openr,unit,file
  endif else begin
     print,'get_log error: source =',source,$
           ' should be a unit number or a filename.'
     retall
  end

  if not keyword_set(verbose) then verbose = 0
; If verbose is a string set the index string to it
  if size(verbose,/type) eq 7 then index=verbose else index=''

; Use buffers for efficient reading
  line  = ''
  nheadline = 0
  isheader  = 1
  headlines = strarr(1)
  buf   = long(10000)
  dbuf  = long(10000)
  nt    = long(0)
  while not eof(unit) do begin
     on_ioerror,close_file

     if isheader then begin
        readf, unit, line

                                ; check if the line contains any character that is not a number
        isheader = 0
        for i = 0, strlen(line)-1 do begin
           if strmatch(strmid(line,i,1), '[!	0123456789dDeE \.+-]') $
           then begin
              isheader = 1
              break
           endif
        endfor

                                ; check if line contains a single number only
        if not isheader then begin
           n = 0
           string_to_array,line, numbers, n
           if n le 1 then isheader=1
        endif
        
        if isheader then begin
           if nheadline eq 0 then $
              headlines(0) = line $
           else $
              headlines = [headlines, line]
           nheadline = nheadline + 1
        endif else begin
                                ; split line into numbers
           string_to_array,line, numbers, nwlog
                                ; create arrays to read data into
           wlog_ = dblarr(nwlog)
           wlog  = dblarr(nwlog,buf)

                                ; read first line
           reads, line, wlog_
           if total(finite(wlog_)) eq nwlog then begin
              wlog(*,0) = wlog_
              nt = 1L
           endif

                                ; find variable names in the header lines
           for i = nheadline - 1, 0, -1 do begin
              line = headlines[i]
              char = strlowcase(strmid(strtrim(line,1),0,1))
              if char ge 'a' and char le 'z' then begin

                                ; Overwrite #START with spaces if present
                 j = strpos(line,'#START')
                 if j ge 0 then strput, line, '      ', j

                                ; split line into names
                 string_to_array, line, wlognames, nname

                                ; if number of names agree we are done
                 if nname eq nwlog then BREAK
              endif
           endfor

           if n_elements(wlognames) ne nwlog then begin
              wlognames = strarr(nwlog)
              for i = 0, nwlog - 1 do $
                 wlognames[i] = 'var'+string(i, format='(i2.2)')
           endif

           if verbose then begin
              if filesource then print,'logfile',index,'  =',file
              print,'headlines',index,':'
              print, format='(a)', strtrim(headlines,2)
              for i=0, nwlog-1 do $
                 print,FORMAT='("  wlog",A,"(*,",I2,")= ",A)',index,i,wlognames(i)
           endif

        endelse
     endif else begin
        readf, unit, wlog_
        if total(finite(wlog_)) eq nwlog then begin
           wlog(*,nt) = wlog_
           nt=nt+1
        endif
        if nt ge buf then begin
           buf=buf+dbuf
           wlog=[[wlog],[dblarr(nwlog,buf)]]
        endif
     endelse

  endwhile
  close_file: if filesource then close,unit

  if verbose then print,'Number of recorded timesteps: nt=',nt

  wlog = transpose(wlog(*,0:nt-1))

  logtime = log_time(wlog,wlognames,timeunit)
  if verbose then print,'Setting logtime',index

end

;=============================================================================

pro plot_log

; Plot variables listed in the space separated func string from the
; files listed in the string array logfilenames.
; Use wlog0...wlognames2 if all the needed arguments are present 
; otherwise read the arrays from the files.
; Use xrange to set the time range, shift times by timeshifts(nlog), 
; set the function ranges with the optional yranges(2,nfunc) array,
; set the colors, the linestyles and the symbols with the 
; colors(nlog), linestyles(nlog) and symbols(nlog) arrays, respectively.
; Smooth data with a boxcar filter of width smooths(nlog). 
; Set the title (default is the list of file names)
; the X title (default is Time [hr]) and 
; the Y titles (defaults are the function names) with the 
; title, xtitle and ytitles(nfunc) arrays, respectively.
; Show the legends (default are the file names) at the position given by
; legendpos array (xmin, xmax, ymin, ymax in a [0,1]x[0,1] box).
; Set the optional variables to zero to get the default behavior.

  common debug_param & on_error, onerror

  common log_data, $
     timeunit, $
     wlog0, logtime0, wlognames0, $ ; renamed from wlog, logtime, wlognames
     wlog1, logtime1, wlognames1, $
     wlog2, logtime2, wlognames2, $
     wlog3, logtime3, wlognames3, $
     wlog4, logtime4, wlognames4, $
     wlog5, logtime5, wlognames5, $
     wlog6, logtime6, wlognames6, $
     wlog7, logtime7, wlognames7, $
     wlog8, logtime8, wlognames8, $
     wlog9, logtime9, wlognames9

  common getlog_param
  common plotlog_param
  common plot_param             ; noerase

  nlog = n_elements(logfilenames)

  string_to_array,logfunc,logfuncs,nlogfunc

; read in arrays if not present
  if    (nlog eq 1 and n_elements(wlognames0) eq 0) $
     or (nlog eq 2 and n_elements(wlognames1) eq 0) $
     or (nlog eq 3 and n_elements(wlognames2) eq 0) then begin

     if nlog ge 1 then get_log,logfilenames(0),wlog0,wlognames0,/verbose
     if nlog ge 2 then get_log,logfilenames(1),wlog1,wlognames1,verbose='1'
     if nlog ge 3 then get_log,logfilenames(2),wlog2,wlognames2,verbose='2'

  endif

; Shift times by 0 unless defined
  if n_elements(timeshifts) lt nlog then timeshifts = fltarr(nlog)

; Shall we do an FFT transform and plot power spectra?
  if n_elements(dofft) eq 0 then dofft=0

; Do not smooth data unless defined
  if n_elements(smooths) lt nlog then smooths = intarr(nlog)

; Calculate the xrange from the data unless defined
  if n_elements(xrange) ne 2 then begin
     DoXrange = 1
     xrange   = [1e30, -1e30]
  endif else $
     DoXrange = 0

; Calculate yranges from the data unless defined
  if n_elements(yranges) ne 2*nlogfunc then begin
     DoYrange = 1
     yranges = fltarr(2,nlogfunc)
     yranges(0,*) =  1e30
     yranges(1,*) = -1e30
  endif else $
     DoYrange = 0

; If line styles are not defined use normal lines (0)
  if n_elements(linestyles) lt nlog then linestyles = intarr(nlog)

; If symbols are not defined do not use symbols (0)
  if n_elements(symbols) lt nlog then symbols = intarr(nlog)

; If colors are not defined, make nlog colors defined to white/black
  if n_elements(colors) lt nlog then colors = intarr(nlog) + 255

; If none of colors, linestyles or symbols are defined, make colors different
  if max(linestyles) eq 0 and max(symbols) eq 0 and min(colors) eq 255 then $
     colors = [255,100,250,150,200,50,25,220,125]

; Define default title
  if n_elements(title) eq 1 and size(title,/type) eq 7 then title0=title $
  else if n_elements(legendpos) ne 4 then title0=strjoin(logfilenames,' ') $
  else                                    title0=' '

; Define default xtitle
  if n_elements(xtitle) eq 1 and size(xtitle,/type) eq 7 then xtitle0=xtitle $
  else if timeunit eq '1' then xtitle0="!5"+'Time' $
  else if dofft then xtitle0="!5"+'Frequency [1/'+timeunit+']' $
  else xtitle0="!5"+'Time ['+timeunit+']'

; Define default ytitles
  if n_elements(ytitles) eq nlogfunc and size(ytitles,/type) eq 7 then $
     ytitles0=ytitles $
  else if dofft then ytitles0 = "!5"+'Power spectrum of '+logfuncs $
  else ytitles0 = "!5" + logfuncs

  if !p.thick eq 0 then begin
     if strpos(!d.name,'X') gt -1 then thick = 1 else thick = 3
  endif else $
     thick = !p.thick

  if strpos(!d.name,'X') gt -1 then loadct,39 else loadct,40

; Set size of plot windows
  spacex = log_spacex*float(!d.x_ch_size)/float(!d.x_size)
  spacey = log_spacey*float(!d.y_ch_size)/float(!d.y_size)
  set_space, nlogfunc, spacex, spacey, sizes, ny = nlogfunc

  if not keyword_set(noerase) then erase

; The first iteration is used to get the X and Y ranges from the
; data. This can be skipped if both ranges are given explicitly.
  if DoXrange or DoYrange then iter0 = 1 else iter0 = 2

  for iter = iter0, 2 do begin
     for ilog = 0, nlog-1 do begin
        case ilog of 
           0: begin
              wlog=wlog0
              wlognames = wlognames0
           end
           1: begin
              wlog=wlog1
              wlognames = wlognames1
           end
           2: begin
              wlog=wlog2
              wlognames = wlognames2
           end
           3: begin
              wlog=wlog3
              wlognames = wlognames3
           end
           4: begin
              wlog=wlog4
              wlognames = wlognames4
           end
           5: begin
              wlog=wlog5
              wlognames = wlognames5
           end
           6: begin
              wlog=wlog6
              wlognames = wlognames6
           end
           7: begin
              wlog=wlog7
              wlognames = wlognames7
           end
           8: begin
              wlog=wlog8
              wlognames = wlognames8
           end
           9: begin
              wlog=wlog9
              wlognames = wlognames9
           end
        endcase
        
        hour = log_time(wlog,wlognames,timeunit) + timeshifts(ilog)

        if(dofft)then begin
           n = n_elements(hour)
           freqmin = 1/(hour(n-1) - hour(0))
           print,'Frequency spacing=',freqmin
           xcoord = findgen(n)*freqmin
        endif else $
           xcoord = hour

        for ifunc = 0, nlogfunc-1 do begin

           field = double(log_func(wlog, wlognames, logfuncs(ifunc), error))

           if error then begin
              if iter eq 1 then $
                 print,"function ",logfuncs(ifunc), $
                       " was not found in wlog",strtrim(string(ilog),2)
           endif else begin

              if(dofft) then field = abs(fft(field))^2

              if(smooths(ilog) gt 1)then field = smooth(field,smooths(ilog))

              if iter eq 1 then begin
                 if DoXrange then begin
                    xrange[0]   = min( [ xrange[0], xcoord ] )
                    xrange[1]   = max( [ xrange[1], xcoord ] )
                 endif else $
                    field = field( where(xcoord ge xrange[0] and $
                                         xcoord le xrange[1]))
                 if DoYrange then begin
                    yranges[0,ifunc] = min( [ yranges[0,ifunc], field ] )
                    yranges[1,ifunc] = max( [ yranges[1,ifunc], field ] )
                 endif
              endif else begin
                 set_position, sizes, 0, ifunc, posm, /rect
                 posm(0) = posm(0) + 0.05
                 if nlogfunc lt 3 then posm(1) = posm(1) + 0.05/nlogfunc
                 title1  = ''
                 xtitle1 = ''
                 xtickname1 = strarr(60)+' '
                 xstyle=5
                 ystyle=5
                 if ilog eq 0 then begin
                    xstyle=1
                    ystyle=1
                    if ifunc eq 0       then title1  = title0
                    if ifunc eq nlogfunc-1 then begin
                       xtitle1    = xtitle0
                       xtickname1 = !x.tickname
                    endif
                 endif
                 plot, xcoord, field, pos = posm, $
                       xrange = xrange, $
                       yrange = yranges(*,ifunc), $
                       ylog = dofft, $
                       xstyle = xstyle, $
                       ystyle = ystyle, $
                       title  = title1, $
                       xtitle = xtitle1, $
                       xtickname = xtickname1, $
                       ytitle = ytitles0(ifunc), $
                       color = colors(ilog), $
                       psym  = symbols(ilog), $
                       linestyle = linestyles(ilog), $
                       thick = thick, $
                       /noerase

                 if ilog eq nlog-1 then oplot,xrange,[0,0],linestyle=2
                 
              endelse
           endelse
        endfor

        if n_elements(legendpos) eq 4 and iter eq 2 then begin
           ;; get vertical position of legend
           ypos=legendpos(3) - (ilog+0.5)/nlog*(legendpos(3)-legendpos(2))
           ;; draw a line (or point) with the appropriate color/linetype/symbol
           plot,legendpos(0:1),[ypos,ypos],    $
                color     = colors(ilog),         $
                psym      = symbols(ilog),        $
                linestyle = linestyles(ilog),     $
                thick     = thick,                $
                /normal, xrange=[0,1], yrange=[0,1], $
                /noerase, xstyle=-1, ystyle=-1
           
                                ; print out legend or logfile name
           if n_elements(legends) eq nlog then legend=legends(ilog) $
           else                                legend=logfilenames(ilog)
           xyouts,legendpos(1), $
                  ypos - float(!d.y_ch_size)/float(!d.y_size)*!p.charsize/2, $
                  '  '+legend
        endif
     endfor
  endfor 

  if DoXrange then xrange =0
  if DoYrange then yranges=0
  if strpos(!d.name,'X') lt 0 then loadct,39

end
;=============================================================================
pro rms_logfiles,logfilename,varname,tmin=tmin,tmax=tmax,verbose=verbose

; Print the rms deviation between two logfiles for variables in varname.
; If varname is not present, show rms for all variables.

  common debug_param & on_error, onerror

  interpol_logfiles,logfilename,var0,var1,varname,time,tmin=tmin,tmax=tmax,$
                    verbose=verbose
  string_to_array,varname,varnames,nvar
  ntime = n_elements(time)

  print,'var rms(A-B) rsm(A) rms(B)'
  for ivar=0,nvar-1 do $
     print,varnames(ivar),sqrt(total((var0(*,ivar)-var1(*,ivar))^2)/ntime), $
           sqrt(total(var0(*,ivar)^2)/ntime), sqrt(total(var1(*,ivar)^2)/ntime)

end
;============================================================================
pro interpol_logfiles,logfilename,var0,var1,varname,time,tmin=tmin,tmax=tmax,$
                      verbose=verbose

; Interpolate variables between two logfiles for variables in varname.
; If varname is not present, interpolate all variables.

  common debug_param & on_error, onerror

  string_to_array,logfilename,logfilenames,nfile

  get_log, logfilenames(0), wlog0, varnames0, verbose=verbose
  get_log, logfilenames(1), wlog1, varnames1, verbose=verbose
  if not keyword_set(varname) then varname = varnames0
  interpol_log,wlog0,wlog1,var0,var1,varname,varnames0,varnames1,$
               time,tmin=tmin,tmax=tmax

end
;============================================================================
pro interpol_log,wlog0,wlog1,var0,var1,varname,varnames0,varnames1,time,$
                 tmin=tmin,tmax=tmax,timeunit=timeunit

; Interpolate the variables listed in varname to the time of wlog0
; between tmin and tmax. 

  common debug_param & on_error, onerror

  string_to_array,varname,varnames,nvar

  if nvar eq 0 then begin
     print,'Usage: interpol_log, wlog0, wlog1, var0, var1, varnames ', $
           '[,varnames0] [,varnames1] [,time] [,tmin=tmin] [,tmax=tmax]'
     retall
  endif

  if n_elements(varnames0) eq 0 then varnames0 = varnames
  if n_elements(varnames1) eq 0 then varnames1 = varnames0

  nvar0 = n_elements(varnames0)
  nvar1 = n_elements(varnames1)

  sizewlog0=size(wlog0)
  sizewlog1=size(wlog1)

  if sizewlog0(0) ne 2 or sizewlog1(0) ne 2 then begin
     print,'wlog0 and wlog1 must be 2D arrays'
     retall
  endif

  if sizewlog0(2) ne nvar0 then begin
     print,'Second dimension of wlog0 should be nvar0=',nvar0
     retall
  endif

  if sizewlog1(2) ne nvar1 then begin
     print,'Second dimension of wlog1 should be nvar1=',nvar1
     retall
  endif

  time0 = log_time(wlog0,varnames0,timeunit)
  time1 = log_time(wlog1,varnames1,timeunit)

  if not keyword_set(tmin) then tmin = max([ min(time0), min(time1) ])
  if not keyword_set(tmax) then tmax = min([ max(time0), max(time1) ])

  index0 = where(time0 ge tmin and time0 le tmax)
  time   = time0(index0)
  ntime  = n_elements(time)

  var0   = fltarr(ntime, nvar)
  var1   = fltarr(ntime, nvar)

  for ivar = 0, nvar-1 do begin

     name = varnames(ivar)

     field0 = log_func(wlog0, varnames0, name, error0)
     field1 = log_func(wlog1, varnames1, name, error1)

     if error0 then print,'Could not find variable ',name,' in wlog0'
     if error1 then print,'Could not find variable ',name,' in wlog1'

     if not error0 and not error1 then begin

        var0(*,ivar) = field0(index0)
        var1(*,ivar) = interpol(field1,time1,time)

     endif
  endfor

end

;==========================================
pro quit
   exit
end
;==========================================

pro set_space, nb, spacex, spacey, sizes, nx = nx, ny = ny

; Determines the size and multiplication factors for plotting perfect circles
; or squares. This routine is used to simply find the parameters, another
; procedure, set_position, is used to actually find the position of the
; circle or square.
; This routine maximizes the area used by the plots, determinining the best
; positions for the number of plots that the user has selected.
;
; input parameters:
; nb - number of plots on a page
; spacex - amount of horizontal space between plots in normalized coordinates
; spacey - amount of vertical   space between plots in normalized coordinates
;
; output parameters:
; bs - box size (size of the plotting region)
; nbx, nby - number of plots in the x and y directions
; xoff, yoff - x and y offsets for positions
; xf, yf - x and y multiplication factors for making perfect squares
;
; This has been adapted to allow the user to define how many objects
;   are in the x and y direction on Jan 2, 1998


  common debug_param & on_error, onerror

  sizes = {bs:0.0, nbx:0, nby:0, xoff:0.0, yoff:0.0, xf:0.0, yf:0.0, $
           ppp: nb, spacex:spacex, spacey:spacey}

  xsi = float(!d.x_size)
  ysi = float(!d.y_size)

  xs = xsi - 5.0*spacex*xsi
  ys = ysi - 5.0*spacey*ysi

  if nb eq 1 then begin

     sizes.nbx = 1
     sizes.nby = 1
     sizes.bs = 1.0 - spacex

     if xs gt ys then begin

        sizes.yf = 1.0
        sizes.xf = ys/xs

     endif else begin

        sizes.xf = 1.0
        sizes.yf = xs/ys

     endelse

  endif else begin

     if (n_elements(nx) gt 0) then begin
        sizes.nbx = nx(0)
        if n_elements(ny) eq 0 then sizes.nby = nb/nx(0) else sizes.nby = ny(0)
     endif else begin
        if (n_elements(ny) gt 0) then begin
           sizes.nby = ny(0)
           sizes.nbx = nb/ny(0)
        endif else begin
           if xs gt ys then begin
              sizes.nbx = round(sqrt(nb))
              sizes.nby = fix(nb/sizes.nbx)
           endif else begin
              sizes.nby = round(sqrt(nb))
              sizes.nbx = fix(nb/sizes.nby)
           endelse
        endelse
     endelse

     if xs gt ys then begin

        if (sizes.nbx*sizes.nby lt nb) then                               $
           if (sizes.nbx le sizes.nby) then sizes.nbx = sizes.nbx + 1      $
           else sizes.nby = sizes.nby + 1                                  $
           else                                                        	$
              if (sizes.nbx lt sizes.nby) and					$
           (n_elements(nx) eq 0) and					$
           (n_elements(ny) eq 0) then begin
           temp = sizes.nby
           sizes.nby = sizes.nbx
           sizes.nbx = temp
        endif

        sizes.yf = 1.0
        sizes.xf = ys/xs
        ;; box size is total size - spacing / number of boxes
        ;; first set this based on X direction
        sizes.bs = ((1.0-spacex*(sizes.nbx-1))/sizes.nbx )/sizes.xf
        ;; check if it fits in the Y direction. If not, set again.
        if sizes.nby*sizes.bs+spacey*(sizes.nby-1) gt 1.0 then 		$
           sizes.bs = (1.0- spacey*(sizes.nby-1))/sizes.nby

     endif else begin

        if (sizes.nbx*sizes.nby lt nb) then				$
           if (sizes.nby le sizes.nbx) then sizes.nby = sizes.nby + 1	$
           else sizes.nbx = sizes.nbx + 1					$
           else								$
              if (sizes.nby lt sizes.nbx) and					$
           (n_elements(nx) eq 0) and					$
           (n_elements(ny) eq 0) then begin
           temp = sizes.nby
           sizes.nby = sizes.nbx
           sizes.nbx = temp
        endif

        sizes.xf = 1.0
        sizes.yf = xs/ys
        sizes.bs = ((1.0 - spacey*(sizes.nby-1))/sizes.nby)/sizes.yf

        if sizes.nbx*sizes.bs+spacex*(sizes.nbx-1) gt 1.0 then $
           sizes.bs = (1.0 - spacex*(sizes.nbx-1))/sizes.nbx

     endelse

  endelse

  sizes.xoff = (1.0 - sizes.xf*(sizes.bs*sizes.nbx + spacex*(sizes.nbx-1)))/2.0
  sizes.yoff = (1.0 - sizes.yf*(sizes.bs*sizes.nby + spacey*(sizes.nby-1)))/2.0

end

;============================================================================

pro set_position, sizes, xipos, yipos, pos, rect = rect, $
                  xmargin = xmargin, ymargin = ymargin

; used in conjunction with set_space. Determines the position of the current
; plotting region, given the output parameters from set_space.
;
; Input parameters:
; nb, spacex, spacey, bs, nbx, nby, xoff, yoff, xf, yf - Outputs from set_space
; pos_num - the number of the plot, ranges from 0 : bs-1
;
; Output parameters:
;
; pos - the position of the plot, used in the plot command
;
; modified to make rectangles on Jan 2, 1998

  common debug_param & on_error, onerror

  nb = sizes.ppp
  spacex = sizes.spacex
  spacey = sizes.spacey

  yf2 = sizes.yf
  yf  = sizes.yf*(1.0-spacey)
  xf2 = sizes.xf
  xf  = sizes.xf*(1.0-spacex)

  if keyword_set(rect) then begin

    if keyword_set(xmargin) then xmar = xmargin(0) 			$
    else xmar = spacex/2.0

    if keyword_set(ymargin) then ymar = ymargin(0) 			$
    else ymar = spacey/2.0

    xbuffer = 3.0*float(!d.x_ch_size)/float(!d.x_size) * !p.charsize +spacex/4.0
    xtotal = 1.0 - (spacex*float(sizes.nbx-1) + xmar + xf2*spacex/2.0) - xbuffer
    xbs = xtotal/(float(sizes.nbx)*xf)

    xoff = xmar - xf2*spacex/2.0 + xbuffer - spacex/4.0

    ybuffer = 3.0*float(!d.y_ch_size)/float(!d.y_size) * !p.charsize
    ytotal = 1.0 - (spacey*float(sizes.nby-1) + ymar + yf2*spacey/2.0) - ybuffer
    ybs = ytotal/(float(sizes.nby)*yf)

    yoff = spacey/4.0

  endif else begin

    xbs  = sizes.bs
    xoff = sizes.xoff
    ybs  = sizes.bs
    yoff = sizes.yoff

  endelse

  ;xpos0 = float(xipos) * (xbs+spacex)*xf + xoff + xf2*spacex/2.0
  ;xpos1 = float(xipos) * (xbs+spacex)*xf + xoff + xf2*spacex/2.0 + xbs*xf

  xpos0 = float(xipos) * (xbs+spacex)*xf + xoff + xf2*spacex
  xpos1 = float(xipos) * (xbs+spacex)*xf + xoff + xf2*spacex + xbs*xf

  ypos0 = (1.0-yf2*spacey/2) - (yipos * (ybs+spacey)*yf + ybs*yf) - yoff
  ypos1 = (1.0-yf2*spacey/2) - (yipos * (ybs+spacey)*yf) - yoff

  pos= [xpos0,ypos0,xpos1,ypos1]

end

;============================================================================
pro plot_color_bar, pos, maxmin

; plot color bar based on the current color table

  common debug_param & on_error, onerror

  xrange=!x.range & yrange=!y.range & !x.range=0 & !y.range=0

  maxi = max(maxmin)
  mini = min(maxmin)

  array = findgen(10,256)
  for i=0,9 do array(i,*) = findgen(256)/(256-1)*(maxi-mini) + mini

  levels=(findgen(60)-1)/(58-1)*(maxi-mini)+mini

;; The !5 in the title makes sure that the fonts produced later
;; will look the same every time when saved into eps/ps file.

  contour, array, /noerase, /cell_fill, xstyle = 5, ystyle = 5, $
           levels = levels, pos=pos, title='!5 '

  plot, maxmin, /noerase, pos = pos, xstyle=1, ystyle=1, /nodata,$
        xtickname = [' ',' '], xticks = 1, xminor=1  , $
        ytickname = strarr(60) + ' ', yticklen = 0.25, $
        title=' ', xtitle=' ',ytitle=' '
  axis, 1, ystyle=1, /nodata, yax=1, charsize=0.9*(!p.charsize > 1.), $
        ytitle=' ', ytickname = strarr(60)

  !x.range=xrange & !y.range=yrange

end

;============================================================================
pro makect, color

  ;; Create color table corresponding to color='mid','blue','red','rwb','bwr'

  common debug_param & on_error, onerror


  common colors

  ;; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

  if not keyword_set(color) then begin

     print,'red   - white to red'
     print,'blue  - white to blue'
     print,'rwb   - red white blue'
     print,'bwr   - blue white red'
     print,'mid   - blue green white yellow red'

     color = ''
     read,'Enter color table from list above : ', color

  endif

  color = strlowcase(color)

  ;; Set read, green, blue to values normalized to the 0.0 -- 1.0 range.

  case color of
     'red' : begin
        r(*) = 1.
        g(*) = 1. - findgen(n)/(n-1)
        b(*) = 1. - findgen(n)/(n-1)
     end

     'blue' : begin
        r(*) = 1. - findgen(n)/(n-1)
        b(*) = 1.
        g(*) = 1. - findgen(n)/(n-1)
     end

     'rwb' : begin
        half=n/2
        r(0:half-1) = 1.
        g(0:half-1) = findgen(half)/(half-1)
        b(0:half-1) = findgen(half)/(half-1)

        r(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        b(half:n-1) = 1.
     end

     'bwr' : begin
        half=n/2
        b(0:half-1) = 1.
        g(0:half-1) = findgen(half)/(half-1)
        r(0:half-1) = findgen(half)/(half-1)
        
        b(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
        r(half:n-1) = 1.
     end

     'mid' : begin
        r(0:n/3-1)     = 0.0
        r(n/3:n/2-1)   = findgen(n/2-n/3)/(n/2-n/3-1)
        r(n/2:n-1)     = 1.0
        
        b(0:n/2-1)      = 1.
        b(n/2:2*n/3-1)  = 1. - findgen(2*n/3-n/2)/(2*n/3-n/2-1)
        b(2*n/3-1:n-1)  = 0.
        
        g(0:n/3-1)      = findgen(n/3)/(n/3-1)
        g(n/3:2*n/3-1)  = 1.
        g(2*n/3:n-1)    = 1. - findgen(n-2*n/3)/(n-2*n/3-1)
        
     end

     else : begin
        print, "Unknown value for color=",color
        r(*) = findgen(n)
        g(*) = findgen(n)
        b(*) = findgen(n)
     end

  endcase

  r(0) = 0.0
  g(0) = 0.0
  b(0) = 0.0

  r(n-1) = 1.0
  g(n-1) = 1.0
  b(n-1) = 1.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

;=============================================================================
pro save_pict, filename, headline, varname, w, x, $
               it, time, eqpar, ndim=ndim, gencoord=gencoord, $
               filetype=filetype, append=append

  common debug_param & on_error, onerror

  if n_elements(filename) eq 0 or n_elements(headline) eq 0 or $
     n_elements(varname) eq 0 or n_elements(w) eq 0 then begin
     print,'ERROR in save_pict: ', $
           'arguments unit headline, varname, w are required'
     retall
  endif

  if n_elements(filetype) eq 0 then filetype = 'ascii'
  if n_elements(it) eq 0 then it=0
  if n_elements(time) eq 0 then time=0.0
  neqpar = n_elements(eqpar)

  sw = size(w)

  if n_elements(x) gt 0 then begin
     sx = size(x)
     if n_elements(ndim) eq 0 then ndim = sx(0) - 1 > 1
     nx = sx(1:ndim)
  endif else begin
     if n_elements(ndim) eq 0 then ndim = sw(0) - 1 > 1
     nx = sw(1:ndim)
     case ndim of
        1: x = findgen(nx(0)) + 1
        2: begin
           x = fltarr(nx(0), nx(1), 2)
           for j = 0L, nx(1)-1 do for i = 0L, nx(0)-1 do x(i,j,*) = [i+1, j+1]
        end
        3: begin
           x = fltarr(nx(0), nx(1), nx(2)) + 1
           for k = 0L, nx(2)-1 do for j = 0L, nx(1)-1 do $
              for i = 0L, nx(0)-1 do x(i,j,*) = [i+1., j+1., k+1.]
        end
     endcase
  endelse

  ; number of variables
  if sw[0] eq ndim + 1 then nw = sw(ndim+1) else nw = 1

  if keyword_set(gencoord) then ndim = -ndim

  unit=1

  if filetype eq 'ascii' then begin
     openw, unit, filename, append=append

     printf, unit, headline
     printf, unit, it, time, ndim, neqpar, nw, format='(i8, 1e13.5, 3i3)'
     printf, unit, nx, format='(3i8)'
     if neqpar gt 0 then printf, unit, eqpar, format='(100(1e13.5))'
     printf, unit, varname
     case abs(ndim) of
        1: for i=0L, nx(0)-1 do $
           printf, unit, x(i), w(i,*), format='(100(1e18.10))'
        2: for j =0L, nx(1)-1 do for i=0L, nx(0)-1 do $
           printf, unit, x(i,j,*), w(i,j,*), format='(100(1e18.10))'
        3: for k=0L, nx(2)-1 do for j=0L, nx(1)-1 do for i=0L, nx(0)-1 do $
           printf, unit, x(i,j,k,*), w(i,j,k,*), format='(100(1e18.10))'
     endcase
  endif else begin
     ; extend strings to 500 characters
     for i = 1, 500-strlen(headline) do headline = headline + ' '
     for i = 1, 500-strlen(varname)  do varname  = varname  + ' '

     ; convert reals to 4 or 8 bytes
     if filetype eq 'real4' then begin
        time  = float(time)
        if neqpar gt 0 then eqpar = float(eqpar)
        x     = float(x)
        w     = float(w)
     endif else begin
        time  = double(time)
        if neqpar gt 0 then eqpar = double(eqpar)
        x     = double(x)
        w     = double(w)
     endelse

     ; write unformatted Fortran file
     openw, unit, filename, /f77_unformatted, append=append
     writeu, unit, headline
     writeu, unit, long(it), time, long(ndim), long(neqpar), long(nw)
     writeu, unit, long(nx)
     if neqpar gt 0 then writeu, unit, eqpar
     writeu, unit, varname
     writeu, unit, x
     case abs(ndim) of
        1: for iw=0, nw-1 do writeu, unit, w(*,iw)
        2: for iw=0, nw-1 do writeu, unit, w(*,*,iw)
        3: for iw=0, nw-1 do writeu, unit, w(*,*,*,iw)
     endcase
  endelse
  close,unit

end

;=============================================================================
pro save_log, filename, headline, varname, array, format=format

  common debug_param & on_error, onerror

  unit=1
  close, unit

  s = size(array)
  ndim = s(0)
  if ndim gt 2 then begin
     print, 'Error in save_log: array should be 1 or 2D array'
     return
  end

  n = s(1)
  print, 'save_log: number of data lines to write is ',n

  if not keyword_set(format) then format = '(100(e13.5))'

  openw, unit, filename
  printf, unit, headline
  printf, unit, varname
  for i = 0, n-1 do printf, unit, array(i,*), format=format
  close, unit

end
;==============================================================================
pro reset_axis

  common debug_param & on_error, onerror

  !x.tickname = strarr(60)
  !y.tickname = strarr(60)

end
