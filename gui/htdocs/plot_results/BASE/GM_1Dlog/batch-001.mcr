#!MC 1000
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

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

### x axis variable
$!LINEMAP [1]  ASSIGN{XAXISVAR = 2}
$!VIEW AXISNICEFIT
  AXIS = 'X' 
  AXISNUM = 1

### y axis variable
$!LINEMAP [1]  ASSIGN{YAXISVAR = 13}
$!VIEW AXISNICEFIT
  AXIS = 'Y' 
  AXISNUM = 1

### set range (X axis)
#$!XYLINEAXIS XDETAIL 1 {RANGEMIN = }
#$!XYLINEAXIS XDETAIL 1 {RANGEMAX = }

### set range (Y axis)
#$!XYLINEAXIS YDETAIL 1 {RANGEMIN = }
#$!XYLINEAXIS YDETAIL 1 {RANGEMAX = }

#$!ATTACHTEXT 
#  XYPOS
#    {
#    X = 20.
#    Y = 84.
#    }
#  TEXTSHAPE
#    {
#    HEIGHT = 24
#    }
#  ATTACHTOZONE = NO
#  ANCHOR = LEFT
#  TEXT = ''

### save file
$!PAPER ORIENTPORTRAIT = YES
$!PRINTSETUP PALETTE = COLOR
$!PRINTSETUP SENDPRINTTOFILE = YES
$!PRINTSETUP PRINTFNAME = 'print.cps'
$!PRINT 
