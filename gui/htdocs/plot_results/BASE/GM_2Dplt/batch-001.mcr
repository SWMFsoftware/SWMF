#!MC 1000

### set useful constants
$!Varset |PI| = (2.*asin(1.))
$!Varset |d2r| = (|PI|/180.)
$!Varset |r2d| = (180./|PI|)

### apply style
$!READSTYLESHEET  "style.sty" 
  INCLUDEPLOTSTYLE = YES
  INCLUDETEXT = YES
  INCLUDEGEOM = YES
  INCLUDEAUXDATA = YES
  INCLUDESTREAMPOSITIONS = YES
  INCLUDECONTOURLEVELS = YES
  MERGE = NO
  INCLUDEFRAMESIZEANDPOSITION = YES

$!RUNMACROFUNCTION  "Set 0,0,0 Rotate Origin" 
$!RUNMACROFUNCTION  "Z=constant" 

### circle at origin
$!ATTACHGEOM 
  FILLCOLOR = BLACK
  ISFILLED = YES
  GEOMTYPE = CIRCLE
  RAWDATA
2.7

### turn on grid on slices
$!FIELD [1]  MESH{COLOR = BLACK}
$!FIELD [1]  MESH{SHOW = YES}
$!FIELD [1]  MESH{LINETHICKNESS = 0.1}

### variable to plot
$!GLOBALCONTOUR 1  VAR = 4

### reset contours
$!RUNMACROFUNCTION  "Reset Contours (MIN/MAX)" 

### set manual contour range
#$!CONTOURLEVELS NEW
#  RAWDATA
#1
#0.
#$!LOOP 200
#  $!VarSet |ContToAdd| = (0. + (|LOOP| * ((1.-0.)/200.) ) )
#  $!CONTOURLEVELS ADD
#    RAWDATA
#  1
#  |ContToAdd|
#$!ENDLOOP

### reposition contour legend
$!GLOBALCONTOUR 1  LEGEND{XYPOS{X = 84.5}}
$!GLOBALCONTOUR 1  LEGEND{XYPOS{Y = 89.5}}

$!ATTACHTEXT 
  XYPOS
    {
    X = 10.
    Y = 88.
    }
  TEXTSHAPE
    {
    HEIGHT = 24
    }
  ATTACHTOZONE = NO
  ANCHOR = LEFT
  TEXT = '&(AUXZONE[1]:TIMESIM)'

### set view area
$!TWODAXIS XDETAIL{RANGEMIN = -32}
$!TWODAXIS XDETAIL{RANGEMAX = 16}
$!TWODAXIS YDETAIL{RANGEMIN = -24}
$!TWODAXIS YDETAIL{RANGEMAX = 24}

### vectors
#$!GLOBALTWODVECTOR UVAR = 8
#$!GLOBALTWODVECTOR VVAR = 9
#$!STREAMATTRIBUTES MAXSTEPS = 3000
#$!RESETVECTORLENGTH 
#$!VarSet |SpacingX| = (5.*abs((16)-(-32))/48.)
#$!VarSet |SpacingY| = (5.*abs((24)-(-24))/48.)
#$!VarSet |Yval| = (|SpacingY| * int(-24 / |SpacingY|))
#$!WHILE |Yval| < 24
#  $!VarSet |Xoffset| = (abs(|Yval| / |SpacingY|))
#  $!WHILE |Xoffset| > |SpacingX|
#    $!Varset |Xoffset| = (|Xoffset| - |SpacingX|)
#  $!ENDWHILE
#  $!VarSet |Xval| = ((|SpacingX| * int((-32+|SpacingY|) / |SpacingX|)) - |Xoffset|)
#  $!WHILE |Xval| < 16
#    $!Varset |DIST| = (abs(|Xval|) + abs(|Yval|))
#    $!IF |DIST| > 6.
#      $!STREAMTRACE ADD
#        STREAMTYPE = TWODLINE
#        DIRECTION = BOTH
#        STARTPOS
#          {
#          X = |Xval|
#          Y = |Yval|
#          }
#    $!ENDIF
#    $!Varset |Xval| += |SpacingX|
#  $!ENDWHILE
#  $!Varset |Yval| += |SpacingY|
#$!ENDWHILE
#$!GLOBALSTREAM COLOR = BLACK

### save file
$!PAPER ORIENTPORTRAIT = YES
$!PRINTSETUP PALETTE = COLOR
$!PRINTSETUP SENDPRINTTOFILE = YES
$!PRINTSETUP PRINTFNAME = 'print.cps'
$!PRINT 
