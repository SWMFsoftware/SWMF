#!MC 1000
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

##################################################################
#                                                                #
#            Default tecplot.mcr file.                           #
#                                                                #
#  This file is processed automatically by tecplot on startup.   #
#  The macro functions defined here will appear in the quick     #
#  macro panel (from the "Tools" menu).                          #
#                                                                #
##################################################################

$!VARSET |debug| = 0
$!Varset |PI| = (2.*asin(1.))
$!Varset |d2r| = (|PI|/180.)
$!Varset |r2d| = (180./|PI|)

########################################################### PAGE 1
########################################################### PAGE 1
########################################################### PAGE 1


$!Macrofunction Name = "Clear all text"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!ATTACHTEXT 
    TEXT = ' ' 
  $!PICK ADDALL
    SELECTTEXT = YES
  $!PICK CLEAR
  $!DRAWGRAPHICS TRUE
  $!REDRAW
$!EndMacroFunction


$!Macrofunction Name = "Clear all geometries"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!ATTACHGEOM 
    COLOR = WHITE
    GEOMTYPE = RECTANGLE
    RAWDATA
  0.001 0.001
  $!PICK ADDALL
    SELECTGEOMS = YES
  $!PICK CLEAR
  $!DRAWGRAPHICS TRUE
  $!REDRAW
$!EndMacroFunction


$!VarSet |LegendOnOff| = 0
$!Macrofunction Name = "Toggle Contour Legend"
  ShowInMacroPanel = True
  $!VarSet |LegendOnOff| += 1
  $!IF |LegendOnOff| >= 2
    $!VarSet |LegendOnOff| = 0
  $!ENDIF
  $!IF |LegendOnOff| == 0
    $!GLOBALCONTOUR LEGEND{SHOW = NO}
  $!ENDIF
  $!IF |LegendOnOff| == 1
    $!RUNMACROFUNCTION "Set Contour Legend" ("20")
  $!ENDIF
  $!REDRAW
$!Endmacrofunction


$!Macrofunction Name = "Set Contour Legend"
  ShowInMacroPanel = False
  $!GLOBALCONTOUR LEGEND{SHOW = YES}
  $!GLOBALCONTOUR LEGEND{OVERLAYBARGRID = NO}
  $!GLOBALCONTOUR LABELS{AUTOLEVELSKIP = |1|}
  $!GLOBALCONTOUR LEGEND{TEXTSHAPE{HEIGHT = 2}}
  $!GLOBALCONTOUR LEGEND{BOX{COLOR = WHITE}}
  $!GLOBALCONTOUR LEGEND{BOX{BOXTYPE = NONE}}
  $!GLOBALCONTOUR LEGEND{LABELLOCATION = COLORMAPDIVISIONS}
$!Endmacrofunction


$!Macrofunction Name = "Add Circle at Origin"
  ShowInMacroPanel = True
  $!PromptForTextString |CircleRadius|
    Instructions = "Enter radius of circle."
  $!ATTACHGEOM 
    FILLCOLOR = BLACK
    ISFILLED = YES
    GEOMTYPE = CIRCLE
    RAWDATA
  |CircleRadius|
$!Endmacrofunction


$!VarSet |ActiveZone| = 0
$!MacroFunction Name = "Set Zone: First"
  ShowInMacroPanel = True
  $!VarSet |ActiveZone| = 1
  $!ActiveFieldZones = [|ActiveZone|]
$!EndMacroFunction


$!MacroFunction Name = "Set Zone: Next"
  ShowInMacroPanel = True
  $!Varset |ActiveZone| = (|ActiveZone| + |ActiveZoneStride|)
  $!if |ActiveZone| > |NumZones|
    $!VarSet |ActiveZone| = |NumZones|
  $!endif
  $!if |ActiveZone| < 1
    $!VarSet |ActiveZone| = 1
  $!endif
  $!ActiveFieldZones = [|ActiveZone|]
$!EndMacroFunction


$!Varset |ActiveZoneStride| = 1
$!MacroFunction Name = "Next / Prev  Stride:"
  ShowInMacroPanel = True
  $!PromptForTextString |ActiveZoneStride|
    Instructions = "Enter the stride for Next/Prev zones [|ActiveZoneStride|]"
$!EndMacroFunction


$!MacroFunction Name = "Set Zone: Prev"
  ShowInMacroPanel = True
  $!Varset |ActiveZone| = (|ActiveZone| - |ActiveZoneStride|)
  $!if |ActiveZone| > |NumZones|
    $!VarSet |ActiveZone| = |NumZones|
  $!endif
  $!if |ActiveZone| < 1
    $!VarSet |ActiveZone| = 1
  $!endif
  $!ActiveFieldZones = [|ActiveZone|]
$!EndMacroFunction


$!MacroFunction Name = "Set Zone: Last"
  ShowInMacroPanel = True
  $!VarSet |ActiveZone| = |NumZones|
  $!ActiveFieldZones = [|ActiveZone|]
$!EndMacroFunction


$!MacroFunction Name = "Label Zones: Number"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!LOOP |NUMZONES|
    $!ATTACHTEXT 
      XYPOS
        {
        X = 1
        Y = 98.5
        }
      ZONE = |Loop|
      ATTACHTOZONE = YES
      ANCHOR = HEADLEFT
      TEXT = 'Zone #|Loop%03d|'
  $!ENDLOOP
  $!DRAWGRAPHICS TRUE
$!EndMacroFunction


$!MacroFunction Name = "Label Zones: Time"
  ShowInMacroPanel = True
  $!PromptForTextString |ZoneFirst|
    Instructions = "Enter first zone number to label.\n -1 for all zones."
  $!If |ZoneFirst| > 0
    $!PromptForTextString |ZoneLast|
      Instructions = "Enter last zone number. (-1 for last zone)"
    $!If |ZoneLast| == -1
      $!Varset |ZoneLast| = |NUMZONES|
    $!Endif
  $!Endif
  $!If |ZoneFirst| == -1
    $!Varset |ZoneFirst| = 1
    $!Varset |ZoneLast| = |NUMZONES|
  $!Endif
  $!PromptForTextString |TimeH|
    Instructions = "Enter start time hour."
  $!PromptForTextString |TimeM|
    Instructions = "Enter start time minute."
  $!PromptForTextString |TimeS|
    Instructions = "Enter start time second."
  $!PromptForTextString |TimeInc|
    Instructions = "Enter time increment in seconds."
  $!PromptForTextString |TimeSuffix|
    Instructions = "Enter time label suffix (space for nothing)."
  $!VarSet |LoopZones| = (|ZoneLast| -|ZoneFirst| + 1)
  $!LOOP |LoopZones|
    $!Varset |CurrentZone| = (|ZoneFirst| + |Loop| - 1)
    $!WHILE |TimeS| >= 60
      $!Varset |TimeS| = (|TimeS| - 60)
      $!Varset |TimeM| = (|TimeM| + 1)
    $!ENDWHILE
    $!WHILE |TimeM| >= 60
      $!Varset |TimeM| = (|TimeM| - 60)
      $!Varset |TimeH| = (|TimeH| + 1)
    $!ENDWHILE
    $!ATTACHTEXT 
      XYPOS
        {
        X = 1
        Y = 1.5
        }
      TEXTSHAPE
        {
        HEIGHT = 28
        }
      ZONE = |CurrentZone|
      ATTACHTOZONE = YES
      ANCHOR = LEFT
      TEXT = 'T = |TimeH%02d|:|TimeM%02d|:|TimeS%02d| |TimeSuffix|'
    $!Varset |TimeS| = (|TimeS| + |TimeInc|)
  $!ENDLOOP
$!EndMacroFunction


$!MacroFunction Name = "Label Zones: AUXZONE"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!PROMPTFORTEXTSTRING |DataName|
    INSTRUCTIONS = "Enter Aux Data zone 'Data Name'."
  $!LOOP |NUMZONES|
    $!ATTACHTEXT 
      XYPOS
        {
        X = 5
        Y = 5
        }
      TEXTSHAPE
        {
        HEIGHT = 20
        }
      ZONE = |Loop|
      ATTACHTOZONE = YES
      ANCHOR = HEADLEFT
      TEXT = '&(AUXZONE[|Loop|]:|DataName|)'
  $!ENDLOOP
  $!DRAWGRAPHICS TRUE
$!EndMacroFunction


#$!MacroFunction Name = "Reset Contours (Nice)"
#  ShowInMacroPanel = True
#  $!DRAWGRAPHICS FALSE
#  $!CONTOURLEVELS RESETTONICE
#    CONTOURGROUP = 1
#    APPROXNUMVALUES = 101
#  $!DRAWGRAPHICS TRUE
#  $!REDRAW 
#$!EndMacroFunction


$!MacroFunction Name = "Reset Contours (MIN/MAX)"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!RUNMACROFUNCTION "Reset Contours" ("201")
  $!DRAWGRAPHICS TRUE
  $!REDRAW 
$!EndMacroFunction


$!VarSet |SavedNCont| = 201
$!MacroFunction Name = "Reset Contours"
  ShowInMacroPanel = False
  $!VarSet |NCont| = |1|
  $!If |NCont| == 0
    $!VarSet |NCont| = |SavedNCont|
  $!Endif
  $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = |MinC|}}
  $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = |MaxC|}}
  $!ContourLevels New
     RawData
      2
       |MinC|
       |MaxC|
  $!Loop |NCont|
     $!VarSet |Coeff2| = ( (|loop|-1)/(|NCont|-1) )
     $!VarSet |Coeff1| = ( (|NCont| - |loop|)/(|NCont|-1) )
     $!VarSet |CurrentLevel| = (|Coeff1|*|MinC| + |Coeff2|*|MaxC|)
     $!ContourLevels Add
        RawData
         1
         |CurrentLevel|
  $!EndLoop
  $!RemoveVar |NCont|
  $!RemoveVar |Coeff1|
  $!RemoveVar |Coeff2|
  $!RemoveVar |CurrentLevel|
  $!VarSet |LegendOnOff| = 1
  $!RUNMACROFUNCTION "Set Contour Legend" ("20")
$!EndMacroFunction


$!VarSet |ContMin| =  -1.0
$!VarSet |ContMax| =   1.0
$!VarSet |ContNum| = 201
$!MacroFunction Name = "Set and Save Contours"
  ShowInMacroPanel = True
  $!PROMPTFORTEXTSTRING |ContMin|
    INSTRUCTIONS = "Enter minimum contour value."
  $!PROMPTFORTEXTSTRING |ContMax|
    INSTRUCTIONS = "Enter maximum contour value."
  $!PROMPTFORTEXTSTRING |ContNum|
    INSTRUCTIONS = "Enter number of contours."
  $!RUNMACROFUNCTION "Reapply Saved Contours"
$!EndMacroFunction


$!MacroFunction Name = "Reapply Saved Contours"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!CONTOURLEVELS NEW
    RAWDATA
  1
  |ContMin|
  $!VarSet |ContLoop| = (|ContNum| - 1)
  $!LOOP |ContLoop|
    $!VarSet |ContToAdd| = (|ContMin| + (|LOOP| * ( (|ContMax|-|ContMin|) / (|ContNum|-1)) ) )
    $!CONTOURLEVELS ADD
      RAWDATA
    1
    |ContToAdd|
  $!ENDLOOP
  $!RUNMACROFUNCTION "Set Contour Legend" ("20")
  $!DRAWGRAPHICS TRUE
  $!REDRAW 
$!EndMacroFunction


########################################################### PAGE 2
########################################################### PAGE 2
########################################################### PAGE 2


$!MacroFunction Name = "CSEM Logo"
  ShowInMacroPanel = True
  $!ATTACHGEOM 
    GEOMTYPE = GEOMIMAGE
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 0.5
      Y = 90
      }
    IMAGEFILENAME = '/usr/local/CSEM/CSEM-circle-trans.png' 
    PIXELASPECTRATIO = 1
    RESIZEFILTER = GAUSSIANFILTER
    RAWDATA
      8.51923078299 9.00000001453 
  $!ATTACHTEXT 
    ANCHORPOS
      {
      X = 10
      Y = 98
      }
    TEXTSHAPE
      {
      HEIGHT = 14
      }
    ANCHOR = HEADLEFT
    TEXT = 'Center for Space Environment Modeling\nUniversity of Michigan' 
  $!REDRAW
$!EndMacroFunction


$!Macrofunction Name = "Save Temporary Style"
  ShowInMacroPanel = True
  $!WRITESTYLESHEET  "/tmp/originaltecplot.sty" 
    INCLUDEPLOTSTYLE = YES
    INCLUDETEXT = YES
    INCLUDEGEOM = YES
    INCLUDEAUXDATA = YES
    INCLUDESTREAMPOSITIONS = YES
    INCLUDECONTOURLEVELS = YES
    INCLUDEFACTORYDEFAULTS = NO
    USERELATIVEPATHS = NO
    COMPRESS = NO
$!Endmacrofunction


$!Macrofunction Name = "Load Temporary Style"
  ShowInMacroPanel = True
  $!READSTYLESHEET  "/tmp/originaltecplot.sty" 
    INCLUDEPLOTSTYLE = YES
    INCLUDETEXT = YES
    INCLUDEGEOM = YES
    INCLUDEAUXDATA = YES
    INCLUDESTREAMPOSITIONS = YES
    INCLUDECONTOURLEVELS = YES
    MERGE = NO
    INCLUDEFRAMESIZEANDPOSITION = NO
$!Endmacrofunction


$!Varset |SaveRestoreStyle| = 1
$!MacroFunction Name = "Toggle Macro S/L Style"
  $!Varset |SaveRestoreStyle| = (abs(|SaveRestoreStyle|-1))
  $!If |SaveRestoreStyle| == 0
    $!Pause "A style file will be saved at beginning of macros and can be recovered."
  $!EndIf
  $!If |SaveRestoreStyle| == 1
    $!Pause "No style file will be saved"
  $!EndIf
$!EndMacroFunction


$!VarSet |XRev| = 0
$!Macrofunction Name = "Reverse X Axis (2D)"
  ShowInMacroPanel = True
  $!VarSet |XRev| += 1
  $!IF |XRev| >= 2
    $!VarSet |XRev| = 0
  $!ENDIF
  $!IF |XRev| == 0
    $!TWODAXIS XDETAIL{ISREVERSED = NO}
  $!ENDIF
  $!IF |XRev| == 1
    $!TWODAXIS XDETAIL{ISREVERSED = YES}
  $!ENDIF
  $!REDRAW
$!Endmacrofunction


$!VarSet |YRev| = 0
$!Macrofunction Name = "Reverse Y Axis (2D)"
  ShowInMacroPanel = True
  $!VarSet |YRev| += 1
  $!IF |YRev| >= 2
    $!VarSet |YRev| = 0
  $!ENDIF
  $!IF |YRev| == 0
    $!TWODAXIS YDETAIL{ISREVERSED = NO}
  $!ENDIF
  $!IF |YRev| == 1
    $!TWODAXIS YDETAIL{ISREVERSED = YES}
  $!ENDIF
  $!REDRAW
$!Endmacrofunction


$!Macrofunction Name = "Set 0,0,0 Rotate Origin"
  ShowInMacroPanel = True
  $!Varset |CenterX| = 0.
  $!Varset |CenterY| = 0.
  $!Varset |CenterZ| = 0.
  $!GLOBALTHREED ROTATEORIGIN{X = |CenterX|}
  $!GLOBALTHREED ROTATEORIGIN{Y = |CenterY|}
  $!GLOBALTHREED ROTATEORIGIN{Z = |CenterZ|}
$!Endmacrofunction


$!MacroFunction Name = "+/- Contours"
  ShowInMacroPanel = True
  $!CONTOURLEVELS NEW
    RAWDATA
  3
  -.001
  0
  .001
  $!RUNMACROFUNCTION "Set Contour Legend" ("1")
  $!REDRAW 
$!EndMacroFunction


$!MacroFunction Name = "Black Grids/Boundaries"
  ShowInMacroPanel = True
  $!FIELD [1-|Numzones|]  MESH{COLOR = BLACK}
  $!FIELD [1-|Numzones|]  BOUNDARY{COLOR = BLACK}
  $!REDRAW
$!EndMacroFunction


$!MacroFunction Name = "Ionosphere 2D"
  ShowInMacroPanel = True
  $!FIELDLAYERS SHOWBOUNDARY = NO
  $!FIELDLAYERS SHOWMESH = NO
  $!FRAMEMODE = TWOD
  $!GLOBALCONTOUR VAR = 8
  $!FIELDLAYERS SHOWCONTOUR = YES
  $!TWODAXIS XVAR = 2
  $!TWODAXIS YVAR = 1
  $!TWODAXIS XDETAIL{ISREVERSED = YES}
  $!TWODAXIS XDETAIL{SHOWAXIS = NO}
  $!TWODAXIS YDETAIL{SHOWAXIS = NO}
  $!TWODAXIS GRIDAREA{EXTENTS{X1 = 0}}
  $!TWODAXIS GRIDAREA{EXTENTS{Y2 = 100}}
  $!TWODAXIS GRIDAREA{EXTENTS{Y1 = 0}}
  $!TWODAXIS GRIDAREA{EXTENTS{X2 = 100}}
  $!TWODAXIS XDETAIL{AXISPOSITION = 0}
  $!TWODAXIS YDETAIL{AXISPOSITION = 0}
  $!TWODAXIS XDETAIL{RANGEMIN = -1.30632911392}
  $!TWODAXIS XDETAIL{RANGEMAX = .8}
  $!TWODAXIS YDETAIL{RANGEMIN = -.8}
  $!TWODAXIS YDETAIL{RANGEMAX = .8}
  $!LOOP |NUMZONES|
    $!ATTACHTEXT 
      POSITIONCOORDSYS = GRID
      XYPOS
        {
        X = -.2
        Y = -.67
        }
      TEXTSHAPE
        {
        SIZEUNITS = GRID
        HEIGHT = 0.07
        }
      ANCHOR = HEADLEFT
      ZONE = |Loop|
      ATTACHTOZONE = YES
      TEXT = '&(ZONENAME:|Loop|)' 
  $!ENDLOOP
  $!ATTACHTEXT 
    POSITIONCOORDSYS = GRID
    XYPOS
      {
      X = 0.
      Y = -.67
      }
    TEXTSHAPE
      {
      SIZEUNITS = GRID
      HEIGHT = 0.07
      }
    ANCHOR = HEADCENTER
    TEXT = '00' 
  $!ATTACHTEXT 
    POSITIONCOORDSYS = GRID
    XYPOS
      {
      X =  .67
      Y = 0.
      }
    TEXTSHAPE
      {
      SIZEUNITS = GRID
      HEIGHT = 0.07
      }
    ANCHOR = MIDRIGHT
    TEXT = '18' 
  $!ATTACHTEXT 
    POSITIONCOORDSYS = GRID
    XYPOS
      {
      X = 0.
      Y =  .67
      }
    TEXTSHAPE
      {
      SIZEUNITS = GRID
      HEIGHT = 0.07
      }
    ANCHOR = CENTER
    TEXT = '12' 
  $!ATTACHTEXT 
    POSITIONCOORDSYS = GRID
    XYPOS
      {
      X = -.67
      Y =  0.
      }
    TEXTSHAPE
      {
      SIZEUNITS = GRID
      HEIGHT = 0.07
      }
    ANCHOR = MIDLEFT
    TEXT = '06' 
  $!ATTACHGEOM 
    GEOMTYPE = CIRCLE
    LINEPATTERN = DASHED
    RAWDATA
  0.173648178 
  $!ATTACHGEOM 
    GEOMTYPE = CIRCLE
    LINEPATTERN = DASHED
    RAWDATA
  0.342020143 
  $!ATTACHGEOM 
    GEOMTYPE = CIRCLE
    LINEPATTERN = DASHED
    RAWDATA
  0.500000000 
  $!ATTACHGEOM 
    GEOMTYPE = CIRCLE
    RAWDATA
  0.642787610 
#  $!GLOBALCONTOUR LEGEND{SHOW = YES}
#  $!GLOBALCONTOUR LEGEND{OVERLAYBARGRID = NO}
#  $!GLOBALCONTOUR LEGEND{BOX{BOXTYPE = NONE}}
#  $!GLOBALCONTOUR LEGEND{TEXTSHAPE{HEIGHT = 3}}
#  $!GLOBALCONTOUR LEGEND{XYPOS{Y = 80}}
#  $!GLOBALCONTOUR LEGEND{XYPOS{X = 82}}
  $!GLOBALCONTOUR 1 LEGEND{SHOW = YES}
  $!GLOBALCONTOUR 1 LEGEND{OVERLAYBARGRID = NO}
  $!GLOBALCONTOUR 1 LEGEND{BOX{BOXTYPE = NONE}}
  $!GLOBALCONTOUR 1 LEGEND{TEXTSHAPE{HEIGHT = 3}}
  $!GLOBALCONTOUR 1 LEGEND{XYPOS{X = 94}}
  $!GLOBALCONTOUR 1 LEGEND{XYPOS{Y = 82.5}}
  $!CONTOURLEVELS NEW
    CONTOURGROUP = 1
    RAWDATA
  15
  -0.21
  -0.18
  -0.15
  -0.12
  -0.09
  -0.06
  -0.03
  0
  0.03
  0.06
  0.09
  0.12
  0.15
  0.18
  0.21
  $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
  $!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = -0.21}}
  $!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 0.21}}
  $!COLORMAP 
    CONTOURCOLORMAP = LGRAINBOW
  $!COLORMAPCONTROL RESETTOFACTORY
  $!COLORMAP 
    LGRAINBOW
      {
      CONTROLPOINT 1
        {
        COLORMAPFRACTION = 0
        LEADRGB
          {
          R = 0
          G = 0
          B = 255
          }
        TRAILRGB
          {
          R = 0
          G = 0
          B = 255
          }
        }
      CONTROLPOINT 2
        {
        COLORMAPFRACTION = 0.25
        LEADRGB
          {
          R = 0
          G = 255
          B = 255
          }
        TRAILRGB
          {
          R = 0
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 3
        {
        COLORMAPFRACTION = 0.45
        LEADRGB
          {
          R = 255
          G = 255
          B = 255
          }
        TRAILRGB
          {
          R = 205
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 4
        {
        COLORMAPFRACTION = 0.5
        LEADRGB
          {
          R = 255
          G = 255
          B = 255
          }
        TRAILRGB
          {
          R = 255
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 5
        {
        COLORMAPFRACTION = 0.55
        LEADRGB
          {
          R = 255
          G = 255
          B = 205
          }
        TRAILRGB
          {
          R = 255
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 6
        {
        COLORMAPFRACTION = 0.75
        LEADRGB
          {
          R = 255
          G = 255
          B = 0
          }
        TRAILRGB
          {
          R = 255
          G = 255
          B = 0
          }
        }
      CONTROLPOINT 7
        {
        COLORMAPFRACTION = 1
        LEADRGB
          {
          R = 255
          G = 0
          B = 0
          }
        TRAILRGB
          {
          R = 255
          G = 0
          B = 0
          }
        }
      }
  $!REDRAW 
$!EndMacroFunction


$!MacroFunction Name = "Ionosphere 3D"
  ShowInMacroPanel = True
  $!FIELDLAYERS SHOWBOUNDARY = NO
  $!FIELDLAYERS SHOWMESH = NO
  $!FIELDLAYERS SHOWBOUNDARY = YES
  $!FIELD [1-|NUMZONES|]  BOUNDARY{COLOR = BLACK}
  $!FIELD [1-|NUMZONES|]  BOUNDARY{LINETHICKNESS = 0.4}
  $!FRAMEMODE = THREED
  $!GLOBALCONTOUR VAR = 8
  $!FIELDLAYERS SHOWCONTOUR = YES
  $!THREEDAXIS XDETAIL{SHOWAXIS = NO}
  $!THREEDAXIS YDETAIL{SHOWAXIS = NO}
  $!THREEDAXIS ZDETAIL{SHOWAXIS = NO}
  $!THREEDAXIS FRAMEAXIS{XYPOS{X = 15}}
  $!GLOBALTHREED ROTATEORIGIN{X = 0}
  $!GLOBALTHREED ROTATEORIGIN{Y = 0}
  $!GLOBALTHREED ROTATEORIGIN{Z = 0}
  $!THREEDVIEW VIEWERPOSITION{X = -9.10969208104}
  $!THREEDVIEW VIEWERPOSITION{Y = 3.12268767626}
  $!THREEDVIEW VIEWERPOSITION{Z = 19.6125443991}
  $!THREEDVIEW PSIANGLE = 26.81
  $!THREEDVIEW THETAANGLE = 110.42
  $!THREEDVIEW ALPHAANGLE = -22.65
  $!VIEW SETMAGNIFICATION
    MAG = 1.3
  $!LOOP |NUMZONES|
  $!ATTACHTEXT 
    XYPOS
      {
      X =  95.
      Y =   5.
      }
    ANCHOR = RIGHT
    ZONE = |Loop|
    ATTACHTOZONE = YES
    TEXT = '&(ZONENAME:|Loop|)' 
  $!ENDLOOP
#  $!GLOBALCONTOUR LEGEND{SHOW = YES}
#  $!GLOBALCONTOUR LEGEND{OVERLAYBARGRID = NO}
#  $!GLOBALCONTOUR LEGEND{BOX{BOXTYPE = NONE}}
#  $!GLOBALCONTOUR LEGEND{TEXTSHAPE{HEIGHT = 3}}
#  $!GLOBALCONTOUR LEGEND{XYPOS{Y = 80}}
#  $!GLOBALCONTOUR LEGEND{XYPOS{X = 82}}
  $!GLOBALCONTOUR 1 LEGEND{SHOW = YES}
  $!GLOBALCONTOUR 1 LEGEND{OVERLAYBARGRID = NO}
  $!GLOBALCONTOUR 1 LEGEND{BOX{BOXTYPE = NONE}}
  $!GLOBALCONTOUR 1 LEGEND{TEXTSHAPE{HEIGHT = 3}}
  $!GLOBALCONTOUR 1 LEGEND{XYPOS{X = 94}}
  $!GLOBALCONTOUR 1 LEGEND{XYPOS{Y = 82.5}}
  $!CONTOURLEVELS NEW
    CONTOURGROUP = 1
    RAWDATA
  15
  -0.21
  -0.18
  -0.15
  -0.12
  -0.09
  -0.06
  -0.03
  0
  0.03
  0.06
  0.09
  0.12
  0.15
  0.18
  0.21
  $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
  $!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = -0.21}}
  $!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 0.21}}
  $!COLORMAP 
    CONTOURCOLORMAP = LGRAINBOW
  $!COLORMAPCONTROL RESETTOFACTORY
  $!COLORMAP 
    LGRAINBOW
      {
      CONTROLPOINT 1
        {
        COLORMAPFRACTION = 0
        LEADRGB
          {
          R = 0
          G = 0
          B = 255
          }
        TRAILRGB
          {
          R = 0
          G = 0
          B = 255
          }
        }
      CONTROLPOINT 2
        {
        COLORMAPFRACTION = 0.25
        LEADRGB
          {
          R = 0
          G = 255
          B = 255
          }
        TRAILRGB
          {
          R = 0
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 3
        {
        COLORMAPFRACTION = 0.45
        LEADRGB
          {
          R = 255
          G = 255
          B = 255
          }
        TRAILRGB
          {
          R = 205
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 4
        {
        COLORMAPFRACTION = 0.5
        LEADRGB
          {
          R = 255
          G = 255
          B = 255
          }
        TRAILRGB
          {
          R = 255
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 5
        {
        COLORMAPFRACTION = 0.55
        LEADRGB
          {
          R = 255
          G = 255
          B = 205
          }
        TRAILRGB
          {
          R = 255
          G = 255
          B = 255
          }
        }
      CONTROLPOINT 6
        {
        COLORMAPFRACTION = 0.75
        LEADRGB
          {
          R = 255
          G = 255
          B = 0
          }
        TRAILRGB
          {
          R = 255
          G = 255
          B = 0
          }
        }
      CONTROLPOINT 7
        {
        COLORMAPFRACTION = 1
        LEADRGB
          {
          R = 255
          G = 0
          B = 0
          }
        TRAILRGB
          {
          R = 255
          G = 0
          B = 0
          }
        }
      }
  $!REDRAW 
$!EndMacroFunction


$!Macrofunction Name = "X=constant"
  ShowInMacroPanel = True
  $!TWODAXIS XVAR = 2
  $!TWODAXIS YVAR = 3
  $!TWODAXIS XDETAIL{ISREVERSED = NO}
  $!TWODAXIS YDETAIL{ISREVERSED = NO}
  $!TWODAXIS GRIDAREA{DRAWBORDER = YES}
  $!GLOBALCONTOUR LEGEND{ANCHORALIGNMENT = TOPLEFT}
  $!GLOBALCONTOUR LEGEND{XYPOS{X = 87.4}}
  $!GLOBALCONTOUR LEGEND{XYPOS{Y = 92.7}}
  $!RUNMACROFUNCTION "Black Grids/Boundaries"
$!Endmacrofunction


$!Macrofunction Name = "Y=constant"
  ShowInMacroPanel = True
  $!TWODAXIS XVAR = 1
  $!TWODAXIS YVAR = 3
  $!TWODAXIS XDETAIL{ISREVERSED = YES}
  $!TWODAXIS YDETAIL{ISREVERSED = NO}
  $!TWODAXIS GRIDAREA{DRAWBORDER = YES}
  $!GLOBALCONTOUR LEGEND{ANCHORALIGNMENT = TOPLEFT}
  $!GLOBALCONTOUR LEGEND{XYPOS{X = 87.4}}
  $!GLOBALCONTOUR LEGEND{XYPOS{Y = 92.7}}
  $!RUNMACROFUNCTION "Black Grids/Boundaries"
$!Endmacrofunction


$!Macrofunction Name = "Z=constant"
  ShowInMacroPanel = True
  $!TWODAXIS XVAR = 1
  $!TWODAXIS YVAR = 2
  $!TWODAXIS YDETAIL{ISREVERSED = YES}
  $!TWODAXIS XDETAIL{ISREVERSED = YES}
  $!TWODAXIS GRIDAREA{DRAWBORDER = YES}
  $!GLOBALCONTOUR LEGEND{ANCHORALIGNMENT = TOPLEFT}
  $!GLOBALCONTOUR LEGEND{XYPOS{X = 87.4}}
  $!GLOBALCONTOUR LEGEND{XYPOS{Y = 92.7}}
  $!RUNMACROFUNCTION "Black Grids/Boundaries"
$!Endmacrofunction


$!Macrofunction Name = "RCM 2d"
  ShowInMacroPanel = True
  $!RUNMACROFUNCTION  "Z=constant" 
  $!FIELDLAYERS SHOWBOUNDARY = NO
  $!BLANKING VALUE{INCLUDE = YES}
  $!BLANKING VALUE{BLANKENTIRECELL = NO}
  $!BLANKING VALUE{CONSTRAINT 1 {INCLUDE = YES}}
  $!BLANKING VALUE{CONSTRAINT 1 {VARA = 3}}
  $!BLANKING VALUE{CONSTRAINT 1 {CONSTRAINTOP2MODE = USEVAR}}
  $!BLANKING VALUE{CONSTRAINT 1 {VARB = 8}}
  $!BLANKING VALUE{CONSTRAINT 1 {SHOW = YES}}
  $!TWODAXIS XDETAIL{RANGEMIN = -15}
  $!TWODAXIS YDETAIL{RANGEMIN = -10}
  $!TWODAXIS YDETAIL{RANGEMAX = 10}
$!Endmacrofunction


# $!MacroFunction Name = "RCM rayValues setup"
#   ShowInMacroPanel = True
#   $!TWODAXIS XVAR = 6
#   $!TWODAXIS YVAR = 7
#   $!TWODAXIS XDETAIL{ISREVERSED = YES}
#   $!TWODAXIS YDETAIL{ISREVERSED = YES}
#   $!FIELD [1-|Numzones|]  MESH{COLOR = BLACK}
#   $!FIELD [1-|Numzones|]  BOUNDARY{COLOR = BLACK}
#   $!Varset |CenterX| = 0.
#   $!Varset |CenterY| = 0.
#   $!Varset |CenterZ| = 0.
#   $!GLOBALTHREED ROTATEORIGIN{X = |CenterX|}
#   $!GLOBALTHREED ROTATEORIGIN{Y = |CenterY|}
#   $!GLOBALTHREED ROTATEORIGIN{Z = |CenterZ|}
#   $!View Datafit
#   $!REDRAWALL
# $!EndMacroFunction


########################################################### PAGE 3
########################################################### PAGE 3
########################################################### PAGE 3


$!Macrofunction Name = "Setup Animation"
  ShowInMacroPanel = False
  $!PromptForTextString |ImageWidth|
    Instructions = "Enter the desired image width."
  $!If |Animation| == 1
    $!Varset |Extension| = "avi"
    $!Varset |format| = "AVI"
    $!EXPORTSETUP EXPORTFORMAT = |format|
    $!EXPORTSETUP IMAGEWIDTH = |ImageWidth|
    $!EXPORTSETUP EXPORTFNAME = "|AnimType|.|Extension|"
    $!EXPORTSTART 
    $!EXPORTNEXTFRAME 
  $!Endif
  $!If |Animation| > 1
    $!If |Animation| == 2
      $!Varset |Extension| = "tif"
      $!Varset |format| = "TIFF"
      $!EXPORTSETUP EXPORTFORMAT = |format|
    $!Endif
    $!If |Animation| == 3
      $!Varset |Extension| = "png"
      $!Varset |format| = "PNG"
      $!EXPORTSETUP EXPORTFORMAT = |format|
      $!PRINTSETUP PALETTE = COLOR
    $!Endif
    $!If |Animation| == 4
      $!Varset |Extension| = "ps"
      $!Varset |format| = "PS"
      $!EXPORTSETUP EXPORTFORMAT = |format|
      $!PRINTSETUP PALETTE = COLOR
    $!Endif
    $!PromptForTextString |ImagePath|
      Instructions = "Enter the path to image directory."
    $!EXPORTSETUP IMAGEWIDTH = |ImageWidth|
    $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/|AnimType|0000.|Extension|"
    $!EXPORT
  $!Endif
$!Endmacrofunction


$!Macrofunction Name = "Animate Multiple"
  ShowInMacroPanel = True
  ###
  ### Set defaults
  $!Varset |DoSlice|=0
  $!Varset |DoZones|=0
  $!Varset |DoTranslate|=0
  $!Varset |DoZoom|=0
  $!Varset |DoRotate|=0
  ###
  ### Start macro
  $!If |SaveRestoreStyle| == 1
    $!RUNMACROFUNCTION "Save Temporary Style"
  $!Endif
  $!PromptForTextString |AnimationType|
    Instructions = "Enter the animation type:\n  Slice=16 + Zones=8 + Translate=4 + Zoom=2 + Rotate=1"
  $!If |AnimationType| >= 16
    $!Varset |DoSlice| = 1
    $!Varset |AnimationType| = (|AnimationType| - 16)
  $!Endif
  $!If |AnimationType| >= 8
    $!Varset |DoZones| = 1
    $!Varset |AnimationType| = (|AnimationType| - 8)
  $!Endif
  $!If |AnimationType| >= 4
    $!Varset |DoTranslate| = 1
    $!Varset |AnimationType| = (|AnimationType| - 4)
  $!Endif
  $!If |AnimationType| >= 2
    $!Varset |DoZoom| = 1
    $!Varset |AnimationType| = (|AnimationType| - 2)
  $!Endif
  $!If |AnimationType| >= 1
    $!Varset |DoRotate| = 1
    $!Varset |AnimationType| = (|AnimationType| - 1)
  $!Endif
  $!If |AnimationType| > 0
    $!Pause "WARNING, AnimationType has extra value."
  $!Endif
  ###
  ### Setup zones
  $!If |DoZones| == 1
    $!PromptForTextString |ZoneFirst|
      Instructions = "Enter first zone number.\nOther zones that are active outside of loop will stay active."
    $!Varset |ZoneNumber| = |ZoneFirst|
    $!PromptForTextString |ZoneLast|
      Instructions = "Enter last zone number. (-1 for last zone)"
    $!If |ZoneLast| == -1
      $!Varset |ZoneLast| = |NUMZONES|
    $!Endif
    $!PromptForTextString |VaryContours|
      Instructions = "Do you want the contour value to vary?\n0=NO  1=YES(banded)  2=YES(continuous)"
    $!If |VaryContours| != 0
      $!PromptForTextString |SavedNCont|
        Instructions = "Enter the number of contours."
    $!Endif
    $!PromptForTextString |BzProbe|
      Instructions = "Do you want Bz plotted? 0=NO  1=YES"
    $!ACTIVEFIELDZONES += [|ZoneNumber|]
    $!If |VaryContours| != 0
      $!RUNMACROFUNCTION "Reset Contours" ("0")
      $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = |MinC|}}
      $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = |MaxC|}}
      $!If |VaryContours| == 1
        $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = BANDED}
      $!Endif
      $!If |VaryContours| == 2
        $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
      $!Endif
    $!Endif
    $!If |BzProbe| == 1
      $!Varset |ProbeZone| = |ZoneNumber|
      $!RUNMACROFUNCTION "Probe Bz"
    $!Endif
    $!VarSet |LoopZones| = (|ZoneLast| -|ZoneFirst|)
  $!Endif
  ###
  ### Read animation steps
  $!If |DoZones| == 0
    $!PromptForTextString |NumSteps|
      Instructions = "Enter the number of steps in animation."
  $!Endif
  $!If |DoZones| == 1
    $!PromptForTextString |NumSteps|
      Instructions = "Enter the number of steps in animation.\n NOTE: should be a multiple of number of zones (|LoopZones|)"
  $!Endif
  ###
  ### Setup translation
  $!If |DoTranslate| == 1
    $!PromptForTextString |TransX|
      Instructions = "Enter the percentage of screen to translate in X."
    $!PromptForTextString |TransY|
      Instructions = "Enter the percentage of screen to translate in Y."
    $!Varset |TransStepX| = (|TransX|/|NumSteps|)
    $!Varset |TransStepY| = (|TransY|/|NumSteps|)
  $!Endif
  ###
  ### Setup zoom
  $!If |DoZoom| == 1
    $!PromptForTextString |MagStart|
      Instructions = "Enter the starting zoom magnification.\n NOTE: Set to current value, cancel if unknown."
    $!PromptForTextString |MagFinal|
      Instructions = "Enter the final zoom magnification.\n NOTE: No SlowMo for zoom."
    $!Varset |MagStep| = ((|MagFinal|-|MagStart|)/|NumSteps|)
    $!Varset |MagFactor| = ((|MagFinal|/|MagStart|)**(1./|NumSteps|))
    $!VIEW SETMAGNIFICATION
      MAG = |MagStart|
    $!Varset |MagCurrent| = |MagStart|
  $!Endif
  ###
  ### Setup rotation
  $!If |DoRotate| == 1
    $!PromptForTextString |RotationAxis|
      Instructions = "Enter axis for rotation (X,Y,Z,PSI,...)."
    $!PromptForTextString |RotationTotal|
      Instructions = "Enter total number of degrees to rotate."
    $!Varset |RotationAngle| = (|RotationTotal|/|NumSteps|)
  $!Endif
  ###
  ### Setup slice
  $!If |DoSlice| == 1
    $!If |DoRotate| == 0
      $!PromptForTextString |RotationAxis|
        Instructions = "Enter axis for rotation (X,Y,Z,PSI,...)."
      $!PromptForTextString |RotationTotal|
        Instructions = "Enter total number of degrees to rotate."
      $!Varset |RotationAngle| = (|RotationTotal|/|NumSteps|)
    $!Endif
    $!If "|RotationAxis|" == "x"
      $!Varset |RotationAxis| = "X"
    $!Endif
    $!If "|RotationAxis|" == "y"
      $!Varset |RotationAxis| = "Y" 
    $!Endif
    $!If "|RotationAxis|" == "z"
      $!Varset |RotationAxis| = "Z"
    $!Endif
    $!PromptForTextString |SliceNx|
      Instructions = "Enter X normal component (origin at 0,0,0)."
    $!PromptForTextString |SliceNy|
      Instructions = "Enter Y normal component (origin at 0,0,0)."
    $!PromptForTextString |SliceNz|
      Instructions = "Enter Z normal component (origin at 0,0,0)."
    $!GLOBALTHREED SLICE{ORIGIN{X = 0.0}}
    $!GLOBALTHREED SLICE{ORIGIN{Y = 0.0}}
    $!GLOBALTHREED SLICE{ORIGIN{Z = 0.0}}
    $!If "|RotationAxis|" == "X"
      $!Varset |SliceAngle| = (asin(abs(|SliceNz|)/sqrt(|SliceNy|**2+|SliceNz|**2)))
    $!Endif
    $!If "|RotationAxis|" == "Y"
      $!Varset |SliceAngle| = (asin(abs(|SliceNx|)/sqrt(|SliceNz|**2+|SliceNx|**2)))
    $!Endif
    $!If "|RotationAxis|" == "Z"
      $!Varset |SliceAngle| = (asin(abs(|SliceNy|)/sqrt(|SliceNx|**2+|SliceNy|**2)))
    $!Endif
    $!Varset |COSA| = (cos(|SliceAngle|))
    $!Varset |SINA| = (-sin(|SliceAngle|))
    $!If "|RotationAxis|" == "X"
      $!GLOBALTHREED SLICE{NORMAL{X = 0}}
      $!GLOBALTHREED SLICE{NORMAL{Y = |COSA|}}
      $!GLOBALTHREED SLICE{NORMAL{Z = |SINA|}}
    $!Endif
    $!If "|RotationAxis|" == "Y"
      $!GLOBALTHREED SLICE{NORMAL{X = |SINA|}}
      $!GLOBALTHREED SLICE{NORMAL{Y = 0}}
      $!GLOBALTHREED SLICE{NORMAL{Z = |COSA|}}
    $!Endif
    $!If "|RotationAxis|" == "Z"
      $!GLOBALTHREED SLICE{NORMAL{X = |COSA|}}
      $!GLOBALTHREED SLICE{NORMAL{Y = |SINA|}}
      $!GLOBALTHREED SLICE{NORMAL{Z = 0}}
    $!Endif
    $!CREATESLICEZONEFROMPLANE 
      SLICESOURCE = VOLUMEZONES
      FORCEEXTRACTIONTOSINGLEZONE = YES
    $!Varset |SliceZone| = |NUMZONES|
    ### Format new slice
    $!ACTIVEFIELDZONES += [|SliceZone|]
    $!FIELD [|SliceZone|]  MESH{SHOW = NO}
    $!FIELD [|SliceZone|]  CONTOUR{FLOODCOLORING = GROUP4}
    $!FIELD [|SliceZone|]  SHADE{SHOW = NO}
    $!FIELD [|SliceZone|]  BOUNDARY{SHOW = NO}
    $!FIELD [|SliceZone|]  SURFACEEFFECTS{USETRANSLUCENCY = NO}
    $!Redraw
    $!PromptForTextString |Question|
      Instructions = "Is this correct?  Cancel if not\nNote: Contour colored by C4"
  $!Endif
  ###
  ### Read more parameters
  $!Varset |SlowMo| = 0
  $!If |DoZones| == 0
    $!PromptForTextString |SlowMo|
      Instructions = "Do slower start and stop?\nEnter 0=NO or EVEN number of steps\n(> than NumSteps)."
  $!Endif
  $!PromptForTextString |Animation|
    Instructions = "Enter 0 for Screen, 1 for AVI file\nFrames: 2 for TIFF, 3 for PNG, 4 for PS\n-1 for quick walkthrough"
  ###
  ###
  $!Redraw
  $!If |Animation| == -1
    $!Pause "Initial view."
  $!Endif
  $!If |Animation| > 0
    $!Varset |AnimType| = "Animation"
    $!RUNMACROFUNCTION "Setup Animation"
  $!Endif
  $!If |SlowMo| != 0
    $!Varset |NumSteps| = (|NumSteps| + |SlowMo|)
  $!Endif
  $!If |DoSlice| == 1
    $!DELETEZONES  [|SliceZone|]
  $!Endif
  ###
  ###
  $!Loop |NumSteps|
    $!DRAWGRAPHICS FALSE
    ###
    ### Advance zone
    $!If |DoZones| == 1
      $!Varset |ZoneTmp| = (|ZoneFirst| + int((|Loop| * |LoopZones|) / |NumSteps| ))
      $!If |ZoneTmp| > |ZoneNumber|
        $!ACTIVEFIELDZONES += [|ZoneTmp|]
        $!ACTIVEFIELDZONES -= [|ZoneNumber|]
        $!Varset |ZoneNumber| = |ZoneTmp|
        $!If |VaryContours| != 0
          $!RUNMACROFUNCTION "Reset Contours" ("0")
          $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = |MinC|}}
          $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = |MaxC|}}
        $!Endif
        $!If |BzProbe| == 1
          $!Varset |ProbeZone| = |ZoneNumber|
          $!RUNMACROFUNCTION "Probe Bz"
        $!Endif
      $!Endif
    $!Endif
    ###
    ### Compute SlowMo StepFactor
    $!If |SlowMo| != 0
      $!Varset |StepFactor| = 1.
      $!If |Loop| <= |SlowMo|
        $!Varset |Step| = (|Loop|)
        $!Varset |StepFactor| = ((sin((|Step|*(|PI|/(|SlowMo|+1)))-(0.5*|PI|))+1.)*0.5)
      $!Endif
      $!If |Loop| > (|NumSteps|-|SlowMo|)
        $!Varset |Step| = ((|NumSteps|-|Loop|)+1)
        $!Varset |StepFactor| = ((sin((|Step|*(|PI|/(|SlowMo|+1)))-(0.5*|PI|))+1.)*0.5)
      $!Endif
    $!Endif
    ###
    ### Advance translation
    $!If |DoTranslate| == 1
      $!Varset |ActualTransStepX| = |TransStepX|
      $!Varset |ActualTransStepY| = |TransStepY|
      $!If |SlowMo| != 0
        $!Varset |ActualTransStepX| = (|TransStepX| * |StepFactor|)
        $!Varset |ActualTransStepY| = (|TransStepY| * |StepFactor|)
      $!Endif
      $!VIEW TRANSLATE
        X = |ActualTransStepX|
        Y = |ActualTransStepY|
    $!Endif
    ###
    ### Advance zoom
    $!If |DoZoom| == 9
      $!Varset |ActualMagStep| = |MagStep|
      $!If |SlowMo| != 0
        $!Varset |ActualMagStep| = (|MagStep| * |StepFactor|)
      $!Endif
      $!Varset |MagCurrent| = (|MagCurrent| + |ActualMagStep|)
      $!VIEW SETMAGNIFICATION
        MAG = |MagCurrent|
    $!Endif
    $!If |DoZoom| == 1
      $!Varset |ActualMagFactor| = |MagFactor|
      $!If |SlowMo| != 0
        $!If |Loop| <= |SlowMo|
          $!Varset |ActualMagFactor| = (|ActualMagFactor|**0.5)
        $!Endif
        $!If |Loop| > (|NumSteps|-|SlowMo|)
          $!Varset |ActualMagFactor| = (|ActualMagFactor|**0.5)
        $!Endif
      $!Endif
      $!Varset |MagCurrent| = (|MagCurrent|*|ActualMagFactor|)
      $!VIEW SETMAGNIFICATION
        MAG = |MagCurrent|
    $!Endif
    ###
    ### Advance rotation
    $!If |DoRotate| == 1
      $!Varset |ActualRotation| = |RotationAngle|
      $!If |SlowMo| != 0
        $!Varset |ActualRotation| = (|RotationAngle| * |StepFactor|)
      $!Endif
      $!ROTATE3DVIEW |RotationAxis|
        ANGLE = |ActualRotation|
        ROTATEORIGINLOCATION = DEFINEDORIGIN
    $!Endif
    ###
    ### Advance slice
    $!If |DoSlice| == 1
      $!Varset |ActualRotation| = |RotationAngle|
      $!If |SlowMo| != 0
        $!Varset |ActualRotation| = (|RotationAngle| * |StepFactor|)
      $!Endif
      $!Varset |SliceAngle| += (|ActualRotation|*|d2r|)
      $!Varset |COSA| = (cos(|SliceAngle|))
      $!Varset |SINA| = (-sin(|SliceAngle|))
      $!If "|RotationAxis|" == "X"
        $!GLOBALTHREED SLICE{NORMAL{X = 0}}
        $!GLOBALTHREED SLICE{NORMAL{Y = |COSA|}}
        $!GLOBALTHREED SLICE{NORMAL{Z = |SINA|}}
      $!Endif
      $!If "|RotationAxis|" == "Y"
        $!GLOBALTHREED SLICE{NORMAL{X = |SINA|}}
        $!GLOBALTHREED SLICE{NORMAL{Y = 0}}
        $!GLOBALTHREED SLICE{NORMAL{Z = |COSA|}}
      $!Endif
      $!If "|RotationAxis|" == "Z"
        $!GLOBALTHREED SLICE{NORMAL{X = |COSA|}}
        $!GLOBALTHREED SLICE{NORMAL{Y = |SINA|}}
        $!GLOBALTHREED SLICE{NORMAL{Z = 0}}
      $!Endif
      $!CREATESLICEZONEFROMPLANE 
        SLICESOURCE = VOLUMEZONES
        FORCEEXTRACTIONTOSINGLEZONE = YES
      $!Varset |SliceZone| = |NUMZONES|
      ### Format new slice
      $!ACTIVEFIELDZONES += [|SliceZone|]
      $!FIELD [|SliceZone|]  MESH{SHOW = NO}
      $!FIELD [|SliceZone|]  CONTOUR{FLOODCOLORING = GROUP4}
      $!FIELD [|SliceZone|]  SHADE{SHOW = NO}
      $!FIELD [|SliceZone|]  BOUNDARY{SHOW = NO}
      $!FIELD [|SliceZone|]  SURFACEEFFECTS{USETRANSLUCENCY = NO}
    $!Endif
    ###
    ### Redraw and save frame
    $!DRAWGRAPHICS TRUE
    $!If |Animation| >= 0
      $!Redraw
    $!Endif
    $!If |Animation| == 1
      $!EXPORTNEXTFRAME 
    $!Endif
    $!If |Animation| > 1
      $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/|AnimType||Loop%04d|.|Extension|"
      $!EXPORT
    $!Endif
    $!If |DoSlice| == 1
      $!DELETEZONES  [|SliceZone|]
    $!Endif
  $!Endloop
  ###
  ### Finish up
  $!If |Animation| == -1
    $!Redraw
    $!Pause "Final view."
  $!Endif
  $!If |Animation| == 1
    $!EXPORTFINISH 
  $!Endif
  $!If |SaveRestoreStyle| == 1
    $!PromptForTextString |ReloadStyle|
      Instructions = "Animation is complete. Do you want to return to original view?\n0=NO  1=YES"
    $!If |ReloadStyle| == 1
      $!RUNMACROFUNCTION "Load Temporary Style"
    $!Endif
  $!Endif
$!Endmacrofunction


$!Macrofunction Name = "Animate Zones"
  ShowInMacroPanel = False
  $!If |SaveRestoreStyle| == 1
    $!RUNMACROFUNCTION "Save Temporary Style"
  $!Endif
  $!PromptForTextString |ZoneFirst|
    Instructions = "Enter first zone number.\nOther zones that are active outside of loop will stay active."
  $!PromptForTextString |ZoneLast|
    Instructions = "Enter last zone number. (-1 for last zone)"
  $!If |ZoneLast| == -1
    $!Varset |ZoneLast| = |NUMZONES|
  $!Endif
  $!PromptForTextString |VaryContours|
    Instructions = "Do you want the contour value to vary?\n0=NO  1=YES(banded)  2=YES(continuous)"
  $!If |VaryContours| != 0
    $!PromptForTextString |SavedNCont|
      Instructions = "Enter the number of contours."
  $!Endif
  $!PromptForTextString |BzProbe|
    Instructions = "Do you want Bz plotted? 0=NO  1=YES"
  $!PromptForTextString |Animation|
    Instructions = "Enter 0 for Screen, 1 for AVI file\nFrames: 2 for TIFF, 3 for PNG, 4 for PS"
  $!ACTIVEFIELDZONES += [|ZoneFirst|]
  $!If |VaryContours| != 0
    $!RUNMACROFUNCTION "Reset Contours" ("0")
    $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = |MinC|}}
    $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = |MaxC|}}
    $!If |VaryContours| == 1
      $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = BANDED}
    $!Endif
    $!If |VaryContours| == 2
      $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
    $!Endif
  $!Endif
  $!If |BzProbe| == 1
    $!Varset |ProbeZone| = |ZoneFirst|
    $!RUNMACROFUNCTION "Probe Bz"
  $!Endif
  $!Redraw
  $!If |Animation| != 0
    $!Varset |AnimType| = "ZoneAnimation"
    $!RUNMACROFUNCTION "Setup Animation"
  $!Endif
  $!VarSet |LoopZones| = (|ZoneLast| -|ZoneFirst|)
  $!Loop |LoopZones|
    $!VarSet |ZoneNumber| = (|ZoneFirst| + |Loop|)
    $!ACTIVEFIELDZONES += [|ZoneNumber|]
    $!VarSet |ZoneNumber| = (|ZoneFirst| + |Loop| - 1)
    $!ACTIVEFIELDZONES -= [|ZoneNumber|]
    $!If |VaryContours| != 0
      $!RUNMACROFUNCTION "Reset Contours" ("0")
      $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = |MinC|}}
      $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = |MaxC|}}
    $!Endif
    $!If |BzProbe| == 1
      $!Varset |ProbeZone| = (|ZoneFirst| + |Loop|)
      $!RUNMACROFUNCTION "Probe Bz"
    $!Endif
    $!Redraw
    $!If |Animation| == 1
      $!EXPORTNEXTFRAME 
    $!Endif
    $!If |Animation| > 1
      $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/ZoneAnimation|Loop%04d|.|Extension|"
      $!EXPORT
    $!Endif
  $!Endloop
  $!If |Animation| == 1
    $!EXPORTFINISH 
  $!Endif
  $!If |SaveRestoreStyle| == 1
    $!PromptForTextString |ReloadStyle|
      Instructions = "Animation is complete. Do you want to return to original view?\n0=NO  1=YES"
    $!If |ReloadStyle| == 1
      $!RUNMACROFUNCTION "Load Temporary Style"
    $!Endif
  $!Endif
$!Endmacrofunction


$!Macrofunction Name = "Animate Transparency"
  ShowInMacroPanel = True
  $!If |SaveRestoreStyle| == 1
    $!RUNMACROFUNCTION "Save Temporary Style"
  $!Endif
  $!PromptForTextString |TransZones|
    Instructions = "Enter the zones to vary transparency.\n(Make sure transparency is already on.)"
  $!PromptForTextString |TransStart|
    Instructions = "Enter the starting zone transparency. (1-99)"
  $!Varset |TransStart| = (min(99,max(1,int(|TransStart|))))
  $!PromptForTextString |TransFinal|
    Instructions = "Enter the final zone transparency. (1-99)"
  $!Varset |TransFinal| = (min(99,max(1,int(|TransFinal|))))
  $!PromptForTextString |NumSteps|
    Instructions = "Enter the number of steps to take."
  $!Varset |TransStep| = (int((|TransFinal|-|TransStart|)/|NumSteps|))
  $!Varset |TransCheck| = (|TransStart| + (|NumSteps|*|TransStep|))
  $!If |TransCheck| != |TransFinal|
    $!Pause "WARNING: Steps and endtime don't match.\nStart=|TransStart| nSteps=|NumSteps| Step=|TransStep|\nEnteredEnd=|TransFinal| ComputedEnd=|TransCheck|"
  $!Endif
  $!PromptForTextString |Animation|
    Instructions = "Enter 0 for Screen, 1 for AVI file\nFrames: 2 for TIFF, 3 for PNG, 4 for PS"
  $!FIELD [|TransZones|]  SURFACEEFFECTS{SURFACETRANSLUCENCY = |TransStart|}
  $!Redraw
  $!If |Animation| != 0
    $!Varset |AnimType| = "Transparency"
    $!RUNMACROFUNCTION "Setup Animation"
  $!Endif
  $!Varset |TransCurrent| = |TransStart|
  $!Loop |NumSteps|
    $!Varset |ActualTransStep| = |TransStep|
    $!Varset |TransCurrent| = (|TransCurrent| + |ActualTransStep|)
    $!FIELD [|TransZones|]  SURFACEEFFECTS{SURFACETRANSLUCENCY = |TransCurrent|}
    $!Redraw
    $!If |Animation| == 1
      $!EXPORTNEXTFRAME 
    $!Endif
    $!If |Animation| > 1
      $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/Transparency|Loop%04d|.|Extension|"
      $!EXPORT
    $!Endif
  $!Endloop
  $!If |Animation| == 1
    $!EXPORTFINISH 
  $!Endif
  $!If |SaveRestoreStyle| == 1
    $!PromptForTextString |ReloadStyle|
      Instructions = "Animation is complete. Do you want to return to original view?\n0=NO  1=YES"
    $!If |ReloadStyle| == 1
      $!RUNMACROFUNCTION "Load Temporary Style"
    $!Endif
  $!Endif
$!Endmacrofunction


$!Macrofunction Name = "Animate Flyby"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!PROMPTFORTEXTSTRING |trackzone|
    INSTRUCTIONS = "Enter the zone to use for flyby tracking."
  $!RUNMACROFUNCTION "Save Temporary Style"
  $!ACTIVEFIELDZONES = [|trackzone|]
  $!VARSET |nPoints| = (|MaxI|-1)
  $!RUNMACROFUNCTION "Load Temporary Style"
  $!PROMPTFORTEXTSTRING |fieldofview|
    INSTRUCTIONS = "Enter the field of view angle for flyby."
  $!PROMPTFORTEXTSTRING |pointing|
    INSTRUCTIONS = "Viewing direction?  0.=Path < BLEND < 1.=Origin"
  $!PROMPTFORTEXTSTRING |satellitezone|
    INSTRUCTIONS = "Move zone as satellite?  0=NO, Otherwise zone number"
  $!If |satellitezone| > 0
    $!PROMPTFORTEXTSTRING |satellitetrackzone|
      INSTRUCTIONS = "Zone to track satellite on?  (must mirror flyby tracking zone)"
  $!Endif
  $!PromptForTextString |Animation|
    Instructions = "Enter 0 for Screen, 1 for AVI file\nFrames: 2 for TIFF, 3 for PNG, 4 for PS"
  $!GetFieldValue |xloc|
    Zone = |trackzone|
    Var = 1
    Index = 1
  $!GetFieldValue |yloc|
    Zone = |trackzone|
    Var = 2
    Index = 1
  $!GetFieldValue |zloc|
    Zone = |trackzone|
    Var = 3
    Index = 1
  $!VARSET |dxloc| = |xloc|
  $!VARSET |dyloc| = |yloc|
  $!VARSET |dzloc| = |zloc|
  $!GetFieldValue |dxloc|
    Zone = |trackzone|
    Var = 1
    Index = 2
  $!GetFieldValue |dyloc|
    Zone = |trackzone|
    Var = 2
    Index = 2
  $!GetFieldValue |dzloc|
    Zone = |trackzone|
    Var = 3
    Index = 2
  $!VARSET |dxloc| = (|xloc| - |dxloc|)
  $!VARSET |dyloc| = (|yloc| - |dyloc|)
  $!VARSET |dzloc| = (|zloc| - |dzloc|)
  $!VARSET |xyzloc| = (sqrt( (|xloc|*|xloc|) + (|yloc|*|yloc|) + (|zloc|*|zloc|) ))
  $!VARSET |dxyzloc| = (sqrt( (|dxloc|*|dxloc|) + (|dyloc|*|dyloc|) + (|dzloc|*|dzloc|) ))
  $!VARSET |dxloc| = (((|dxloc|/|dxyzloc|)*(1.-|pointing|)) + ((|xloc|/|xyzloc|)*(|pointing|)))
  $!VARSET |dyloc| = (((|dyloc|/|dxyzloc|)*(1.-|pointing|)) + ((|yloc|/|xyzloc|)*(|pointing|)))
  $!VARSET |dzloc| = (((|dzloc|/|dxyzloc|)*(1.-|pointing|)) + ((|zloc|/|xyzloc|)*(|pointing|)))
  $!DRAWGRAPHICS TRUE
  $!THREEDAXIS XDETAIL{SHOWAXIS = NO}
  $!THREEDAXIS YDETAIL{SHOWAXIS = NO}
  $!THREEDAXIS ZDETAIL{SHOWAXIS = NO}
  $!THREEDVIEW DRAWINPERSPECTIVE = YES
  $!THREEDVIEW FIELDOFVIEW = |fieldofview|
  $!RUNMACROFUNCTION "Set Viewing Perspective"
  $!If |satellitezone| > 0
    $!GetFieldValue |xsat|
      Zone = |satellitetrackzone|
      Var = 1
      Index = 1
    $!GetFieldValue |ysat|
      Zone = |satellitetrackzone|
      Var = 2
      Index = 1
    $!GetFieldValue |zsat|
      Zone = |satellitetrackzone|
      Var = 3
      Index = 1
    $!RUNMACROFUNCTION "Save Temporary Style"
    $!ACTIVEFIELDZONES = [|satellitezone|]
    $!VARSET |satellitePoint| = (|MaxI|)
    $!RUNMACROFUNCTION "Load Temporary Style"
    $!GetFieldValue |xzone|
      Zone = |satellitezone|
      Var = 1
      Index = |satellitePoint|
    $!GetFieldValue |yzone|
      Zone = |satellitezone|
      Var = 2
      Index = |satellitePoint|
    $!GetFieldValue |zzone|
      Zone = |satellitezone|
      Var = 3
      Index = |satellitePoint|
    $!ALTERDATA  [|satellitezone|]
      EQUATION = 'V1=V1+|xsat|-|xzone|' 
    $!ALTERDATA  [|satellitezone|]
      EQUATION = 'V2=V2+|ysat|-|yzone|' 
    $!ALTERDATA  [|satellitezone|]
      EQUATION = 'V3=V3+|zsat|-|zzone|' 
  $!Endif
  $!Redraw
  $!If |Animation| != 0
    $!Varset |AnimType| = "Flyby"
    $!RUNMACROFUNCTION "Setup Animation"
  $!Endif
  $!LOOP |nPoints|
    $!VARSET |looppoint| = (|loop|+1)
    $!GetFieldValue |xloc|
      Zone = |trackzone|
      Var = 1
      Index = |looppoint|
    $!GetFieldValue |yloc|
      Zone = |trackzone|
      Var = 2
      Index = |looppoint|
    $!GetFieldValue |zloc|
      Zone = |trackzone|
      Var = 3
      Index = |looppoint|
    $!VARSET |dxloc| = |xloc|
    $!VARSET |dyloc| = |yloc|
    $!VARSET |dzloc| = |zloc|
    $!IF |looppoint| > |nPoints|
      $!VARSET |looppoint2| = (|looppoint|-1)
      $!GetFieldValue |dxloc|
        Zone = |trackzone|
        Var = 1
        Index = |looppoint2|
      $!GetFieldValue |dyloc|
        Zone = |trackzone|
        Var = 2
        Index = |looppoint2|
      $!GetFieldValue |dzloc|
        Zone = |trackzone|
        Var = 3
        Index = |looppoint2|
      $!VARSET |dxloc| = (|dxloc| - |xloc|)
      $!VARSET |dyloc| = (|dyloc| - |yloc|)
      $!VARSET |dzloc| = (|dzloc| - |zloc|)
    $!ENDIF
    $!IF |looppoint| <= |nPoints|
      $!VARSET |looppoint2| = (|looppoint|+1)
      $!GetFieldValue |dxloc|
        Zone = |trackzone|
        Var = 1
        Index = |looppoint2|
      $!GetFieldValue |dyloc|
        Zone = |trackzone|
        Var = 2
        Index = |looppoint2|
      $!GetFieldValue |dzloc|
        Zone = |trackzone|
        Var = 3
        Index = |looppoint2|
      $!VARSET |dxloc| = (|xloc| - |dxloc|)
      $!VARSET |dyloc| = (|yloc| - |dyloc|)
      $!VARSET |dzloc| = (|zloc| - |dzloc|)
    $!ENDIF
    $!VARSET |xyzloc| = (sqrt( (|xloc|*|xloc|) + (|yloc|*|yloc|) + (|zloc|*|zloc|) ))
    $!VARSET |dxyzloc| = (sqrt( (|dxloc|*|dxloc|) + (|dyloc|*|dyloc|) + (|dzloc|*|dzloc|) ))
    $!VARSET |dxloc| = (((|dxloc|/|dxyzloc|)*(1.-|pointing|)) + ((|xloc|/|xyzloc|)*(|pointing|)))
    $!VARSET |dyloc| = (((|dyloc|/|dxyzloc|)*(1.-|pointing|)) + ((|yloc|/|xyzloc|)*(|pointing|)))
    $!VARSET |dzloc| = (((|dzloc|/|dxyzloc|)*(1.-|pointing|)) + ((|zloc|/|xyzloc|)*(|pointing|)))
    $!RUNMACROFUNCTION "Set Viewing Perspective"
    $!If |satellitezone| > 0
      $!GetFieldValue |xsat|
        Zone = |satellitetrackzone|
        Var = 1
        Index = |looppoint|
      $!GetFieldValue |ysat|
        Zone = |satellitetrackzone|
        Var = 2
        Index = |looppoint|
      $!GetFieldValue |zsat|
        Zone = |satellitetrackzone|
        Var = 3
        Index = |looppoint|
      $!GetFieldValue |xzone|
        Zone = |satellitezone|
        Var = 1
        Index = |satellitePoint|
      $!GetFieldValue |yzone|
        Zone = |satellitezone|
        Var = 2
        Index = |satellitePoint|
      $!GetFieldValue |zzone|
        Zone = |satellitezone|
        Var = 3
        Index = |satellitePoint|
      $!ALTERDATA  [|satellitezone|]
        EQUATION = 'V1=V1+|xsat|-|xzone|' 
      $!ALTERDATA  [|satellitezone|]
        EQUATION = 'V2=V2+|ysat|-|yzone|' 
      $!ALTERDATA  [|satellitezone|]
        EQUATION = 'V3=V3+|zsat|-|zzone|' 
    $!Endif
    $!Redraw
    $!If |Animation| == 1
      $!EXPORTNEXTFRAME 
    $!Endif
    $!If |Animation| > 1
      $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/|AnimType||Loop%04d|.|Extension|"
      $!EXPORT
    $!Endif
  $!Endloop
  $!If |Animation| == 1
    $!EXPORTFINISH
  $!Endif
  $!If |SaveRestoreStyle| == 1
    $!PromptForTextString |ReloadStyle|
      Instructions = "Animation is complete. Do you want to return to original view?\n0=NO  1=YES"
    $!If |ReloadStyle| == 1
      $!RUNMACROFUNCTION "Load Temporary Style"
    $!Endif
  $!Endif
$!Endmacrofunction


$!Varset |DoLighting| = 0
$!Macrofunction Name = "Test viewing perspective"
  ShowInMacroPanel = True
  $!If |SaveRestoreStyle| == 1
    $!RUNMACROFUNCTION "Save Temporary Style"
  $!Endif
  $!PROMPTFORTEXTSTRING |trackzone|
    INSTRUCTIONS = "Enter the zone to use for flyby tracking. (0 to enter point)"
  $!PROMPTFORTEXTSTRING |fieldofview|
    INSTRUCTIONS = "Enter the field of view angle for flyby."
  $!PROMPTFORTEXTSTRING |pointing|
    INSTRUCTIONS = "Viewing direction?  0.=Path < BLEND < 1.=Origin"
  $!PROMPTFORTEXTSTRING |satellitezone|
    INSTRUCTIONS = "Move zone as satellite?  0=NO, Otherwise zone number"
  $!If |satellitezone| > 0
    $!Varset |satellitetrackzone| = 0
    $!If |trackzone| > 0
      $!PROMPTFORTEXTSTRING |satellitetrackzone|
        INSTRUCTIONS = "Zone to track satellite on?  (must mirror flyby tracking zone, 0 to enter point)"
    $!Endif
  $!Endif
  $!If |trackzone| > 0
    $!RUNMACROFUNCTION "Save Temporary Style"
    $!ACTIVEFIELDZONES = [|trackzone|]
    $!VARSET |nPoints| = (|MaxI|-1)
    $!RUNMACROFUNCTION "Load Temporary Style"
    $!GetFieldValue |xloc|
      Zone = |trackzone|
      Var = 1
      Index = 1
    $!GetFieldValue |yloc|
      Zone = |trackzone|
      Var = 2
      Index = 1
    $!GetFieldValue |zloc|
      Zone = |trackzone|
      Var = 3
      Index = 1
    $!VARSET |dxloc| = |xloc|
    $!VARSET |dyloc| = |yloc|
    $!VARSET |dzloc| = |zloc|
    $!GetFieldValue |dxloc|
      Zone = |trackzone|
      Var = 1
      Index = 2
    $!GetFieldValue |dyloc|
      Zone = |trackzone|
      Var = 2
      Index = 2
    $!GetFieldValue |dzloc|
      Zone = |trackzone|
      Var = 3
      Index = 2
    $!VARSET |dxloc| = (|xloc| - |dxloc|)
    $!VARSET |dyloc| = (|yloc| - |dyloc|)
    $!VARSET |dzloc| = (|zloc| - |dzloc|)
    $!VARSET |xyzloc| = (sqrt( (|xloc|*|xloc|) + (|yloc|*|yloc|) + (|zloc|*|zloc|) ))
    $!VARSET |dxyzloc| = (sqrt( (|dxloc|*|dxloc|) + (|dyloc|*|dyloc|) + (|dzloc|*|dzloc|) ))
    $!VARSET |dxloc| = (((|dxloc|/|dxyzloc|)*(1.-|pointing|)) + ((|xloc|/|xyzloc|)*(|pointing|)))
    $!VARSET |dyloc| = (((|dyloc|/|dxyzloc|)*(1.-|pointing|)) + ((|yloc|/|xyzloc|)*(|pointing|)))
    $!VARSET |dzloc| = (((|dzloc|/|dxyzloc|)*(1.-|pointing|)) + ((|zloc|/|xyzloc|)*(|pointing|)))
  $!Endif
  $!If |trackzone| == 0
    $!PROMPTFORTEXTSTRING |xloc|
      INSTRUCTIONS = "Enter the x location for point."
    $!PROMPTFORTEXTSTRING |yloc|
      INSTRUCTIONS = "Enter the y location for point."
    $!PROMPTFORTEXTSTRING |zloc|
      INSTRUCTIONS = "Enter the z location for point."
    $!VARSET |dxloc| = |xloc|
    $!VARSET |dyloc| = |yloc|
    $!VARSET |dzloc| = |zloc|
    $!IF |pointing| < 1.
      $!PROMPTFORTEXTSTRING |dxloc|
        INSTRUCTIONS = "Enter the x location for point 2."
      $!PROMPTFORTEXTSTRING |dyloc|
        INSTRUCTIONS = "Enter the y location for point 2."
      $!PROMPTFORTEXTSTRING |dzloc|
        INSTRUCTIONS = "Enter the z location for point 2."
      $!VARSET |dxloc| = (|xloc| - |dxloc|)
      $!VARSET |dyloc| = (|yloc| - |dyloc|)
      $!VARSET |dzloc| = (|zloc| - |dzloc|)
      $!VARSET |xyzloc| = (sqrt( (|xloc|*|xloc|) + (|yloc|*|yloc|) + (|zloc|*|zloc|) ))
      $!VARSET |dxyzloc| = (sqrt( (|dxloc|*|dxloc|) + (|dyloc|*|dyloc|) + (|dzloc|*|dzloc|) ))
      $!VARSET |dxloc| = (((|dxloc|/|dxyzloc|)*(1.-|pointing|)) + ((|xloc|/|xyzloc|)*(|pointing|)))
      $!VARSET |dyloc| = (((|dyloc|/|dxyzloc|)*(1.-|pointing|)) + ((|yloc|/|xyzloc|)*(|pointing|)))
      $!VARSET |dzloc| = (((|dzloc|/|dxyzloc|)*(1.-|pointing|)) + ((|zloc|/|xyzloc|)*(|pointing|)))
    $!ENDIF
  $!ENDIF
  $!THREEDAXIS XDETAIL{SHOWAXIS = NO}
  $!THREEDAXIS YDETAIL{SHOWAXIS = NO}
  $!THREEDAXIS ZDETAIL{SHOWAXIS = NO}
  $!THREEDVIEW DRAWINPERSPECTIVE = YES
  $!THREEDVIEW FIELDOFVIEW = |fieldofview|
  $!VARSET |debug| = 1
  $!RUNMACROFUNCTION "Set Viewing Perspective"
  $!VARSET |debug| = 0
  $!If |satellitezone| > 0
    $!If |satellitetrackzone| > 0
      $!GetFieldValue |xsat|
        Zone = |satellitetrackzone|
        Var = 1
        Index = 1
      $!GetFieldValue |ysat|
        Zone = |satellitetrackzone|
        Var = 2
        Index = 1
      $!GetFieldValue |zsat|
        Zone = |satellitetrackzone|
        Var = 3
        Index = 1
    $!Endif
    $!If |satellitetrackzone| == 0
      $!PROMPTFORTEXTSTRING |xsat|
        INSTRUCTIONS = "Enter the x location for zone."
      $!PROMPTFORTEXTSTRING |ysat|
        INSTRUCTIONS = "Enter the y location for zone."
      $!PROMPTFORTEXTSTRING |zsat|
        INSTRUCTIONS = "Enter the z location for zone."
    $!ENDIF
    $!RUNMACROFUNCTION "Save Temporary Style"
    $!ACTIVEFIELDZONES = [|satellitezone|]
    $!VARSET |satellitePoint| = (|MaxI|)
    $!RUNMACROFUNCTION "Load Temporary Style"
    $!GetFieldValue |xzone|
      Zone = |satellitezone|
      Var = 1
      Index = |satellitePoint|
    $!GetFieldValue |yzone|
      Zone = |satellitezone|
      Var = 2
      Index = |satellitePoint|
    $!GetFieldValue |zzone|
      Zone = |satellitezone|
      Var = 3
      Index = |satellitePoint|
    $!ALTERDATA  [|satellitezone|]
      EQUATION = 'V1=V1+|xsat|-|xzone|' 
    $!ALTERDATA  [|satellitezone|]
      EQUATION = 'V2=V2+|ysat|-|yzone|' 
    $!ALTERDATA  [|satellitezone|]
      EQUATION = 'V3=V3+|zsat|-|zzone|' 
  $!Endif
  $!Redraw
  $!If |SaveRestoreStyle| == 1
    $!PromptForTextString |ReloadStyle|
      Instructions = "Animation is complete. Do you want to return to original view?\n0=NO  1=YES"
    $!If |ReloadStyle| == 1
      $!RUNMACROFUNCTION "Load Temporary Style"
    $!Endif
  $!Endif
$!Endmacrofunction


### This macro assumes that these variables are defined
###   |xloc| |yloc| and |zloc|
###   |dxloc| |dyloc| and |dzloc|
$!Macrofunction Name = "Set Viewing Perspective"
  ShowInMacroPanel = False
  $!VARSET |dxyloc|  = (sqrt( (|dxloc|*|dxloc|) + (|dyloc|*|dyloc|) ))
  $!VARSET |dxyzloc| = (sqrt( (|dxloc|*|dxloc|) + (|dyloc|*|dyloc|) + (|dzloc|*|dzloc|) ))
  $!VARSET |pdeg| =  ( |r2d| * asin(|dxyloc|/|dxyzloc|) )
  $!IF |dzloc| < 0.
    $!VARSET |pdeg| = (180. - |pdeg|)
  $!ENDIF
  $!IF |dyloc| > 0
    $!VARSET |tdeg| = (( |r2d| * atan(|dxloc|/|dyloc|) ) + 180.)
  $!ENDIF
  $!IF |dyloc| < 0
    $!VARSET |tdeg| = (( |r2d| * atan(|dxloc|/|dyloc|) ) + 360.)
  $!ENDIF
  $!IF |dyloc| == 0
    $!IF |dxloc| >= 0
      $!VARSET |tdeg| = 270.
    $!ENDIF
    $!IF |dxloc| < 0
      $!VARSET |tdeg| = 90.
    $!ENDIF
  $!ENDIF
  $!VARSET |xyzloc| = (sqrt( (|xloc|*|xloc|) + (|yloc|*|yloc|) + (|zloc|*|zloc|) ))
  $!SET3DEYEDISTANCE 
    EYEDISTANCE = |xyzloc|
  $!THREEDVIEW VIEWERPOSITION{X = |xloc|}
  $!THREEDVIEW VIEWERPOSITION{Y = |yloc|}
  $!THREEDVIEW VIEWERPOSITION{Z = |zloc|}
  $!THREEDVIEW ALPHAANGLE = 0
  $!THREEDVIEW PSIANGLE = |pdeg|
  $!THREEDVIEW THETAANGLE = |tdeg|
  $!IF |debug| == 1
    $!PAUSE "xyz: |xloc| |yloc| |zloc|   |xyzloc|\ndxyz: |dxloc| |dyloc| |dzloc|   |dxyloc| |dxyzloc|\nangles: |pdeg| |tdeg|"
  $!ENDIF
  $!RUNMACROFUNCTION "Set Lighting"
$!Endmacrofunction


### This macro assumes that these variables are defined
###   |xloc| |yloc| and |zloc|
$!Macrofunction Name = "Set Lighting"
  ShowInMacroPanel = False
  $!VARSET |xyloc| = (sqrt( (|xloc|*|xloc|) + (|yloc|*|yloc|) ))
  $!VARSET |xyzloc| = (sqrt( (|xloc|*|xloc|) + (|yloc|*|yloc|) + (|zloc|*|zloc|) ))
  $!VARSET |pdeg| =  ( |r2d| * asin(|xyloc|/|xyzloc|) )
  $!IF |zloc| < 0.
    $!VARSET |pdeg| = (180. - |pdeg|)
  $!ENDIF
  $!IF |yloc| > 0
    $!VARSET |tdeg| = (( |r2d| * atan(|xloc|/|yloc|) ) + 180.)
  $!ENDIF
  $!IF |yloc| < 0
    $!VARSET |tdeg| = (( |r2d| * atan(|xloc|/|yloc|) ) + 360.)
  $!ENDIF
  $!IF |yloc| == 0
    $!IF |xloc| >= 0
      $!VARSET |tdeg| = 270.
    $!ENDIF
    $!IF |xloc| < 0
      $!VARSET |tdeg| = 90.
    $!ENDIF
  $!ENDIF
  $!VARSET |tdeg| = (|tdeg| - 180.)
  $!Varset |xlight| = (cos(|d2r|*|tdeg|))
  $!Varset |zlight| = (sin(|d2r|*|tdeg|)*sin(|d2r|*|pdeg|))
  $!Varset |ylight| = (sin(|d2r|*|tdeg|)*cos(|d2r|*|pdeg|))
  $!GLOBALTHREED LIGHTSOURCE{INTENSITY = 75}
  $!GLOBALTHREED LIGHTSOURCE{SURFACECOLORCONTRAST = 100}
  $!GLOBALTHREED LIGHTSOURCE{BACKGROUNDLIGHT = 50}
  $!GLOBALTHREED LIGHTSOURCE{SPECULARSHININESS = 35}
  $!GLOBALTHREED LIGHTSOURCE{SPECULARINTENSITY = 35}
  $!GLOBALTHREED 
    LIGHTSOURCE
      {
      XYZDIRECTION
        {
        X = |xlight|
        Y = |ylight|
        Z = |zlight|
        }
      }
  $!IF |debug| == 1
    $!PAUSE "lighting: |tdeg| |pdeg|\n    |xlight| |ylight| |zlight|"
  $!EndIf
$!Endmacrofunction


$!Varset |ManualFramesNumber| = 0
$!Macrofunction Name = "Export Manual Frames"
  ShowInMacroPanel = True
  $!PromptForTextString |Animation|
    Instructions = "Enter 1 for AVI file\n  AVI file gives warning and popup to do next frame or finish.\nFrames: 2 for TIFF, 3 for PNG, 4 for PS\n  use 'Next Manual Frame' for each desired frame."
  $!Redraw
  $!If |Animation| != 0
    $!Varset |AnimType| = "Frame"
    $!RUNMACROFUNCTION "Setup Animation"
  $!Endif
$!Endmacrofunction


$!Macrofunction Name = "Next Manual Frame"
  ShowInMacroPanel = True
  $!Varset |ManualFramesNumber| += 1
  $!If |Animation| == 1
    $!PAUSE "This macro does nothing for manual AVI frame advance, use popup window."
  $!Endif
  $!If |Animation| > 1
    $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/Frames|ManualFramesNumber%04d|.|Extension|"
    $!EXPORT
  $!Endif
$!Endmacrofunction


$!MacroFunction Name = "Difference two zones"
  ShowInMacroPanel = True
  $!PROMPTFORTEXTSTRING |z1|
    INSTRUCTIONS = "Computing difference of two zones Z1-Z2\nEnter Z1"
  $!PROMPTFORTEXTSTRING |z2|
    INSTRUCTIONS = "Computing difference of two zones Z1-Z2\nEnter Z2"
  $!DUPLICATEZONE 
    SOURCEZONE = |z1|
  $!LINEARINTERPOLATE 
    SOURCEZONES =  [|z2|]
    DESTINATIONZONE = |NUMZONES|
    VARLIST =  [4-|NUMVARS|]
    LINEARINTERPCONST = 0
    LINEARINTERPMODE = SETTOCONST
  $!VarSet |varsloop| = |NUMVARS|
  $!VarSet |varsloop| -= 3
  $!LOOP |varsloop|
    $!VarSet |var| = |Loop|
    $!VarSet |var| += 3
    $!ALTERDATA  [|NUMZONES|]
      EQUATION = 'V|var|=V|var|[|z1|]-V|var|[|NUMZONES|]'
  $!ENDLOOP
  $!ActiveFieldZones = [|NUMZONES|]
  $!FIELD [|NUMZONES|]  MESH{COLOR = BLACK}
  $!FIELD [|NUMZONES|]  BOUNDARY{COLOR = BLACK}
  $!ATTACHTEXT
    XYPOS
      {
      X = 15
      Y = 98.5
      }
    ZONE = |NUMZONES|
    ATTACHTOZONE = YES
    ANCHOR = HEADLEFT
    TEXT = 'Zone Difference: zone |z1| - zone |z2|'
  $!REDRAWALL
$!EndMacroFunction


$!MacroFunction Name = "Probe Bz Upstream"
  ShowInMacroPanel = True
  $!PROMPTFORTEXTSTRING |ProbeZone|
    INSTRUCTIONS = "Enter zone to probe for Bz value.\n NOTE: probes at x=30. y=z=0."
  $!RUNMACROFUNCTION "Probe Bz"
$!EndMacroFunction


$!Varset |ProbeZone| = 1
$!MacroFunction Name = "Probe Bz"
  ShowInMacroPanel = False
  $!DRAWGRAPHICS FALSE
  ###Put something these so that the delete will work
  $!ATTACHGEOM 
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 93.
      Y = 30.
      }
    MACROFUNCTIONCOMMAND = '' 
    RAWDATA
  1
  2
  0 0 
  3. 0.
  ###Delete all text and geometries in that region
  $!PICK ADDALLINRECT
    SELECTGEOMS = YES
    SELECTTEXT = YES
    X1 = 9.
    X2 = 11.
    Y1 = 2.
    Y2 = 8.5
  $!PICK CLEAR
  ###Create zone to interpolate to and get value
  $!CREATERECTANGULARZONE 
    IMAX = 1
    JMAX = 1
    KMAX = 1
    X1 = 30
    Y1 = 0
    Z1 = 0
    X2 = 30
    Y2 = 0
    Z2 = 0
    XVAR = 1
    YVAR = 2
  $!Varset |tmpzone| = |NUMZONES|
  $!ALTERDATA  [|tmpzone|]
    EQUATION = 'V1=30.' 
  $!ALTERDATA  [|tmpzone|]
    EQUATION = 'V2=0.' 
  $!ALTERDATA  [|tmpzone|]
    EQUATION = 'V3=0.' 
  $!INVERSEDISTINTERPOLATE 
    SOURCEZONES =  [|ProbeZone|]
    DESTINATIONZONE = |tmpzone|
    VARLIST =  [4-|NUMVARS|]
    INVDISTEXPONENT = 2.
    INVDISTMINRADIUS = 0
    INTERPPTSELECTION = NEARESTNPOINTS
    INTERPNPOINTS = 4
  $!GetFieldValue |upstreambz|
    Zone = |tmpzone|
    Var = 10
    Index = 1
  $!DELETEZONES  [|tmpzone|]
  ###Draw the lines, text, and vector
  $!ATTACHGEOM 
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 97.
      Y = 30.
      }
    MACROFUNCTIONCOMMAND = '' 
    ARROWHEADSTYLE = FILLED
    ARROWHEADATTACHMENT = ATEND
    ARROWHEADSIZE = 2
    RAWDATA
  1
  2
  0 0 
  0. |upstreambz|
  $!ATTACHGEOM 
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 93.
      Y = 20.
      }
    MACROFUNCTIONCOMMAND = '' 
    RAWDATA
  1
  2
  0 0 
  3. 0.
  $!ATTACHGEOM 
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 93.
      Y = 25.
      }
    MACROFUNCTIONCOMMAND = '' 
    RAWDATA
  1
  2
  0 0 
  3. 0.
  $!ATTACHGEOM 
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 93.
      Y = 30.
      }
    MACROFUNCTIONCOMMAND = '' 
    RAWDATA
  1
  2
  0 0 
  3. 0.
  $!ATTACHGEOM 
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 93.
      Y = 35.
      }
    MACROFUNCTIONCOMMAND = '' 
    RAWDATA
  1
  2
  0 0 
  3. 0.
  $!ATTACHGEOM 
    POSITIONCOORDSYS = FRAME
    ANCHORPOS
      {
      X = 93.
      Y = 40.
      }
    MACROFUNCTIONCOMMAND = '' 
    RAWDATA
  1
  2
  0 0 
  3. 0.
  $!ATTACHTEXT 
    ANCHORPOS
      {
      X = 93
      Y = 30
      }
    ANCHOR = MIDRIGHT
    TEXT = 'Bz=0' 
  $!ATTACHTEXT 
    ANCHORPOS
      {
      X = 93
      Y = 40
      }
    ANCHOR = MIDRIGHT
    TEXT = '10' 
  $!ATTACHTEXT 
    ANCHORPOS
      {
      X = 93
      Y = 35
      }
    ANCHOR = MIDRIGHT
    TEXT = '5' 
  $!ATTACHTEXT 
    ANCHORPOS
      {
      X = 93
      Y = 25
      }
    ANCHOR = MIDRIGHT
    TEXT = '-5' 
  $!ATTACHTEXT 
    ANCHORPOS
      {
      X = 93
      Y = 20
      }
    ANCHOR = MIDRIGHT
    TEXT = '-10' 
  $!DRAWGRAPHICS TRUE
$!EndMacroFunction


$!MacroFunction Name = "Probe point in zones"
  $!PAUSE "This macro will probe a range of zones at a given point\nto a new 1D zone of values to see how a point\nevolves in time."
  $!DRAWGRAPHICS FALSE
  $!PromptForTextString |ZoneFirst|
    Instructions = "Enter first zone number.\nOther zones that are active outside of loop will stay active."
  $!PromptForTextString |ZoneLast|
    Instructions = "Enter last zone number. (-1 for last zone)"
  $!If |ZoneLast| == -1
    $!Varset |ZoneLast| = |NUMZONES|
  $!Endif
  $!VarSet |ZoneFirst| = (|ZoneFirst| - 1)
  $!VarSet |LoopZones| = (|ZoneLast| - |ZoneFirst|)
  $!CREATERECTANGULARZONE 
    IMAX = |LoopZones|
    JMAX = 1
    KMAX = 1
    X1 = 1
    Y1 = 0
    Z1 = 0
    X2 = |LoopZones|
    Y2 = 0
    Z2 = 0
    XVAR = 1
    YVAR = 2
  $!Varset |fillzone| = |NUMZONES|
  $!PromptForTextString |xloc|
    Instructions = "Enter X probe location."
  $!PromptForTextString |yloc|
    Instructions = "Enter Y probe location."
  $!PromptForTextString |zloc|
    Instructions = "Enter Z probe location."
  $!Loop |LoopZones|
    $!Varset |fillloop| = |Loop|
    $!VarSet |ZoneNumber| = (|ZoneFirst| + |Loop|)
    $!ACTIVEFIELDZONES += [|ZoneNumber|]
    $!If |Loop| != 1
      $!VarSet |ZoneNumber| = (|ZoneFirst| + |Loop| - 1)
      $!ACTIVEFIELDZONES -= [|ZoneNumber|]
    $!Endif
    $!Varset |ProbeZone| = (|ZoneFirst| + |Loop|)
    $!CREATERECTANGULARZONE 
      IMAX = 1
      JMAX = 1
      KMAX = 1
      X1 = |xloc|
      Y1 = |yloc|
      Z1 = |zloc|
      X2 = |xloc|
      Y2 = |yloc|
      Z2 = |zloc|
      XVAR = 1
      YVAR = 2
    $!Varset |tmpzone| = |NUMZONES|
    $!INVERSEDISTINTERPOLATE 
      SOURCEZONES =  [|ProbeZone|]
      DESTINATIONZONE = |tmpzone|
      VARLIST =  [4-|NUMVARS|]
      INVDISTEXPONENT = 2.
      INVDISTMINRADIUS = 0
      INTERPPTSELECTION = NEARESTNPOINTS
      INTERPNPOINTS = 4
    $!Loop |NUMVARS|
      $!If |Loop| > 3
        $!GetFieldValue |value|
          Zone = |tmpzone|
          Var = |Loop|
          Index = 1
        $!SetFieldValue
          Zone = |fillzone|
          Var = |Loop|
          Index = |fillloop|
          FieldValue = |value|
      $!Endif
    $!Endloop
    $!DELETEZONES  [|tmpzone|]
  $!Endloop
  $!DRAWGRAPHICS TRUE
$!EndMacroFunction


$!MacroFunction Name = "3D last closed lines"
  $!DRAWGRAPHICS FALSE
  $!PAUSE "Macro assumes:\n  zone 1 is 3D zone\n  variable 15 is 'Status'"
  $!Varset |zone3D| = 1
  $!Varset |StatusVar| = 15
  $!Varset |StatusValue| = 2.9
  $!PromptForTextString |StatusValue|
    Instructions = "Enter the status value [less than 3.].\n (Add 10. to value to change status variable also [result > 5].)"
  $!If |StatusValue| >= 5.
    $!PromptForTextString |StatusVar|
      Instructions = "Enter the new 'Status' variable number"
    $!Varset |StatusValue| = (|StatusValue| - 10.)
  $!Endif
  $!PromptForTextString |ClosedLines|
    Instructions = "Enter the number of last closed lines desired."
  $!PromptForTextString |zloc|
    Instructions = "Enter the z value for search [normally 0., up to +-2]."
  $!GLOBALTHREEDVECTOR UVAR = 8
  $!GLOBALTHREEDVECTOR VVAR = 9
  $!GLOBALTHREEDVECTOR WVAR = 10
  $!RESETVECTORLENGTH 
  $!GLOBALSTREAM COLOR = WHITE
  $!GLOBALSTREAM LINETHICKNESS = 0.4
  $!Loop |ClosedLines|
    ###Find last closed line for a given angle
    $!Varset |Angle| = (((|Loop|-1.)/|ClosedLines|)*(2.*|PI|))
    $!Varset |Radius| = 8.
    $!Varset |delRadius| = 8.
    $!Varset |dir| = 1
    $!Varset |sinA| = (sin(|Angle|))
    $!Varset |cosA| = (cos(|Angle|))
    $!Varset |NumLoops| = 30
    $!Loop |NumLoops|
      $!Varset |xloc| = (|Radius|*|sinA|)
      $!Varset |yloc| = (|Radius|*|cosA|)
      ###Create zone to get status value
      $!CREATERECTANGULARZONE 
        IMAX = 1
        JMAX = 1
        KMAX = 1
        X1 = |xloc|
        Y1 = |yloc|
        Z1 = |zloc|
        X2 = |xloc|
        Y2 = |yloc|
        Z2 = |zloc|
        XVAR = 1
        YVAR = 2
      $!Varset |tmpzone| = |NUMZONES|
      $!INVERSEDISTINTERPOLATE 
        SOURCEZONES =  [|zone3D|]
        DESTINATIONZONE = |tmpzone|
        VARLIST =  [4-|NUMVARS|]
        INVDISTEXPONENT = 2.
        INVDISTMINRADIUS = 0
        INTERPPTSELECTION = NEARESTNPOINTS
        INTERPNPOINTS = 8
      $!GetFieldValue |value|
        Zone = |tmpzone|
        Var = |StatusVar|
        Index = 1
      $!DELETEZONES  [|tmpzone|]
      $!If |value| >= |StatusValue|
        $!If |dir| == -1
          $!Varset |delRadius| = (|delRadius|/(-3.))
	  $!Varset |dir| = 1
        $!Endif
      $!Endif
      $!If |value| < |StatusValue|
        $!If |dir| == 1
          $!Varset |delRadius| = (|delRadius|/(-3.))
	  $!Varset |dir| = -1
        $!Endif
      $!Endif
      $!Varset |Radius| = (|Radius| + |delRadius|)
    $!Endloop
    $!STREAMTRACE ADD
      STREAMTYPE = VOLUMELINE
      DIRECTION = BOTH
      STARTPOS
        {
        X = |xloc|
        Y = |yloc|
        Z = |zloc|
        }
  $!Endloop
  $!DRAWGRAPHICS TRUE
$!EndMacroFunction


$!MacroFunction Name = "CUSTOM Animation"
  ShowInMacroPanel = True
  $!RUNMACROFUNCTION "CUSTOM"
$!EndMacroFunction


$!Varset |Counter| = 0
$!Macrofunction Name = "CUSTOM"
  ShowInMacroPanel = False
  ###
  ### Start macro
  $!If |SaveRestoreStyle| == 1
    $!RUNMACROFUNCTION "Save Temporary Style"
  $!Endif
  ### Loop for CUSTOM steps
  $!Varset |Counter| = 0
  $!Varset |CustomLoop|=0
  $!Loop 5
    $!Varset |CustomLoop| += 1
    ###
    ### Set defaults
    $!Varset |DoZones|=0
    $!Varset |DoTranslate|=0
    $!Varset |DoZoom|=0
    $!Varset |DoRotate|=0
    ###
    $!Varset |AnimationType| = 0
    $!Varset |ZoneFirst| = 0
    $!Varset |ZoneLast| = 0
    $!Varset |VaryContours| = 0
    $!Varset |SavedNCont| = 0
    $!Varset |BzProbe| = 0
    $!Varset |NumSteps| = 0
    $!Varset |TransX| = 0
    $!Varset |TransY| = 0
    $!Varset |MagStart| = 0
    $!Varset |MagFinal| = 0
    $!Varset |RotationAxis| = 0
    $!Varset |RotationTotal| = 0
    $!Varset |SlowMo| = 0
    ###
    $!Varset |SaveFirst| = 0
    $!Varset |DoExtract| = 0
    ###
    $!If |CustomLoop| == 1
      $!OPENLAYOUT  "/Data/ZoomMovie/all2black.lay"
      ###
      $!Varset |AnimationType| = 14
      $!Varset |ZoneFirst| = 1
      $!Varset |ZoneLast| = 33
      $!Varset |VaryContours| = 0
      $!Varset |BzProbe| = 0
      $!Varset |NumSteps| = 32
      $!Varset |TransX| = -30
      $!Varset |TransY| = 0
      $!Varset |MagStart| = 30
      $!Varset |MagFinal| = 15
      ###
      $!Varset |SaveFirst| = 1
      $!Varset |DoExtract| = 1
    $!Endif
    $!If |CustomLoop| == 2
      $!Varset |AnimationType| = 14
      $!Varset |ZoneFirst| = 33
      $!Varset |ZoneLast| = 103
      $!Varset |VaryContours| = 0
      $!Varset |BzProbe| = 0
      $!Varset |NumSteps| = 70
      $!Varset |TransX| = -60
      $!Varset |TransY| = 0
      $!Varset |MagStart| = 15
      $!Varset |MagFinal| = 3
      ###
      $!Varset |SaveFirst| = 0
      $!Varset |DoExtract| = 1
    $!Endif
    $!If |CustomLoop| == 3
      $!Varset |AnimationType| = 14
      $!Varset |ZoneFirst| = 103
      $!Varset |ZoneLast| = 196
      $!Varset |VaryContours| = 0
      $!Varset |BzProbe| = 0
      $!Varset |NumSteps| = 93
      $!Varset |TransX| = -15
      $!Varset |TransY| = 0
      $!Varset |MagStart| = 3
      $!Varset |MagFinal| = 2
      ###
      $!Varset |SaveFirst| = 0
      $!Varset |DoExtract| = 1
    $!Endif
    $!If |CustomLoop| == 4
      $!FIELDLAYERS SHOWMESH = YES
      ###
      $!Varset |AnimationType| = 6
      $!Varset |NumSteps| = 100
      $!Varset |TransX| = -180.4
      $!Varset |TransY| = -48.6
      $!Varset |MagStart| = 2
      $!Varset |MagFinal| = 750
      $!Varset |SlowMo| = 0
      ###
      $!Varset |SaveFirst| = 1
      $!Varset |DoExtract| = 0
    $!Endif
    $!If |CustomLoop| == 5
      $!OPENLAYOUT  "/Data/ZoomMovie/allE2black.lay"
      ###
      $!Varset |AnimationType| = 8
      $!Varset |ZoneFirst| = 1
      $!Varset |ZoneLast| = -1
      $!Varset |VaryContours| = 0
      $!Varset |BzProbe| = 0
      $!Varset |NumSteps| = 110
      ###
      $!Varset |SaveFirst| = 1
      $!Varset |DoExtract| = 1
    $!Endif
#    $!PromptForTextString |AnimationType|
#      Instructions = "Enter the animation type:\n  Zones=8 + Translate=4 + Zoom=2 + Rotate=1"
    $!If |AnimationType| >= 8
      $!Varset |DoZones| = 1
      $!Varset |AnimationType| = (|AnimationType| - 8)
    $!Endif
    $!If |AnimationType| >= 4
      $!Varset |DoTranslate| = 1
      $!Varset |AnimationType| = (|AnimationType| - 4)
    $!Endif
    $!If |AnimationType| >= 2
      $!Varset |DoZoom| = 1
      $!Varset |AnimationType| = (|AnimationType| - 2)
    $!Endif
    $!If |AnimationType| >= 1
      $!Varset |DoRotate| = 1
      $!Varset |AnimationType| = (|AnimationType| - 1)
    $!Endif
    ###
    ### Setup zones
    $!If |DoZones| == 1
#      $!PromptForTextString |ZoneFirst|
#        Instructions = "Enter first zone number.\nOther zones that are active outside of loop will stay active."
      $!Varset |ZoneNumber| = |ZoneFirst|
#      $!PromptForTextString |ZoneLast|
#        Instructions = "Enter last zone number. (-1 for last zone)"
      $!If |ZoneLast| == -1
        $!Varset |ZoneLast| = |NUMZONES|
      $!Endif
#      $!PromptForTextString |VaryContours|
#        Instructions = "Do you want the contour value to vary?\n0=NO  1=YES(banded)  2=YES(continuous)"
      $!If |VaryContours| != 0
#        $!PromptForTextString |SavedNCont|
#          Instructions = "Enter the number of contours."
      $!Endif
#      $!PromptForTextString |BzProbe|
#        Instructions = "Do you want Bz plotted? 0=NO  1=YES"
      $!ACTIVEFIELDZONES += [|ZoneNumber|]
      $!If |VaryContours| != 0
        $!RUNMACROFUNCTION "Reset Contours" ("0")
        $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = |MinC|}}
        $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = |MaxC|}}
        $!If |VaryContours| == 1
          $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = BANDED}
        $!Endif
        $!If |VaryContours| == 2
          $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
        $!Endif
      $!Endif
      $!If |BzProbe| == 1
        $!Varset |ProbeZone| = |ZoneNumber|
        $!RUNMACROFUNCTION "Probe Bz"
      $!Endif
      $!VarSet |LoopZones| = (|ZoneLast| -|ZoneFirst|)
    $!Endif
    ###
    ### Read animation steps
    $!If |DoZones| == 0
#      $!PromptForTextString |NumSteps|
#        Instructions = "Enter the number of steps in animation."
    $!Endif
    $!If |DoZones| == 1
#      $!PromptForTextString |NumSteps|
#        Instructions = "Enter the number of steps in animation.\n NOTE: should be a multiple of number of zones (|LoopZones|)"
    $!Endif
    ###
    ### Setup translation
    $!If |DoTranslate| == 1
#      $!PromptForTextString |TransX|
#        Instructions = "Enter the percentage of screen to translate in X."
#      $!PromptForTextString |TransY|
#        Instructions = "Enter the percentage of screen to translate in Y."
      $!Varset |TransStepX| = (|TransX|/|NumSteps|)
      $!Varset |TransStepY| = (|TransY|/|NumSteps|)
    $!Endif
    ###
    ### Setup zoom
    $!If |DoZoom| == 1
#      $!PromptForTextString |MagStart|
#        Instructions = "Enter the starting zoom magnification.\n NOTE: Set to current value, cancel if unknown."
#      $!PromptForTextString |MagFinal|
#        Instructions = "Enter the final zoom magnification.\n NOTE: No SlowMo for zoom."
      $!Varset |MagStep| = ((|MagFinal|-|MagStart|)/|NumSteps|)
      $!Varset |MagFactor| = ((|MagFinal|/|MagStart|)**(1./|NumSteps|))
      $!VIEW SETMAGNIFICATION
        MAG = |MagStart|
      $!Varset |MagCurrent| = |MagStart|
    $!Endif
    ###
    ### Setup rotation
    $!If |DoRotate| == 1
#      $!PromptForTextString |RotationAxis|
#        Instructions = "Enter axis for rotation (X,Y,Z,PSI,...)."
#      $!PromptForTextString |RotationTotal|
#        Instructions = "Enter total number of degrees to rotate."
      $!Varset |RotationAngle| = (|RotationTotal|/|NumSteps|)
    $!Endif
    ###
    ### Read more parameters
    $!Varset |SlowMo| = 0
    $!If |DoZones| == 0
#      $!PromptForTextString |SlowMo|
#        Instructions = "Do slower start and stop?\nEnter 0=NO or EVEN number of steps\n(> than NumSteps)."
    $!Endif
    ###
    $!Redraw
    ###
    $!If |CustomLoop| == 1
      $!PromptForTextString |Animation|
        Instructions = "Enter 0 for Screen, 1 for AVI file\nFrames: 2 for TIFF, 3 for PNG, 4 for PS\n-1 for quick walkthrough"
      $!If |Animation| == -1
        $!Pause "Initial view."
      $!Endif
      $!If |Animation| > 0
        $!Varset |AnimType| = "Animation"
        $!RUNMACROFUNCTION "Setup Animation"
      $!Endif
      $!If |DoExtract| == 1
        $!RUNMACROFUNCTION "Extract Lines"
      $!Endif
    $!Endif
    $!If |CustomLoop| > 1
      $!If |SaveFirst| == 1
        $!Varset |Counter| += 1
        $!If |Animation| == 1
          $!EXPORTNEXTFRAME 
        $!Endif
        $!If |Animation| > 1
          $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/|AnimType||Counter%04d|.|Extension|"
          $!EXPORT
        $!Endif
      $!Endif
    $!Endif
    ###
    $!If |SlowMo| != 0
      $!Varset |NumSteps| = (|NumSteps| + |SlowMo|)
    $!Endif
    ###
    ###
    $!Loop |NumSteps|
      $!DRAWGRAPHICS FALSE
      ###
      ### Advance zone
      $!If |DoZones| == 1
        $!Varset |ZoneTmp| = (|ZoneFirst| + int((|Loop| * |LoopZones|) / |NumSteps| ))
        $!If |ZoneTmp| > |ZoneNumber|
          $!ACTIVEFIELDZONES += [|ZoneTmp|]
          $!ACTIVEFIELDZONES -= [|ZoneNumber|]
          $!Varset |ZoneNumber| = |ZoneTmp|
          $!If |VaryContours| != 0
            $!RUNMACROFUNCTION "Reset Contours" ("0")
            $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = |MinC|}}
            $!GLOBALCONTOUR COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = |MaxC|}}
          $!Endif
          $!If |BzProbe| == 1
            $!Varset |ProbeZone| = |ZoneNumber|
            $!RUNMACROFUNCTION "Probe Bz"
          $!Endif
        $!Endif
      $!Endif
      ###
      ### Compute SlowMo StepFactor
      $!If |SlowMo| != 0
        $!Varset |StepFactor| = 1.
        $!If |Loop| <= |SlowMo|
          $!Varset |Step| = (|Loop|)
          $!Varset |StepFactor| = ((sin((|Step|*(|PI|/(|SlowMo|+1)))-(0.5*|PI|))+1.)*0.5)
        $!Endif
        $!If |Loop| > (|NumSteps|-|SlowMo|)
          $!Varset |Step| = ((|NumSteps|-|Loop|)+1)
          $!Varset |StepFactor| = ((sin((|Step|*(|PI|/(|SlowMo|+1)))-(0.5*|PI|))+1.)*0.5)
        $!Endif
      $!Endif
      ###
      ### Advance translation
      $!If |DoTranslate| == 1
        $!Varset |ActualTransStepX| = |TransStepX|
        $!Varset |ActualTransStepY| = |TransStepY|
        $!If |SlowMo| != 0
          $!Varset |ActualTransStepX| = (|TransStepX| * |StepFactor|)
          $!Varset |ActualTransStepY| = (|TransStepY| * |StepFactor|)
        $!Endif
        $!VIEW TRANSLATE
          X = |ActualTransStepX|
          Y = |ActualTransStepY|
      $!Endif
      ###
      ### Advance zoom
      $!If |DoZoom| == 9
        $!Varset |ActualMagStep| = |MagStep|
        $!If |SlowMo| != 0
          $!Varset |ActualMagStep| = (|MagStep| * |StepFactor|)
        $!Endif
        $!Varset |MagCurrent| = (|MagCurrent| + |ActualMagStep|)
        $!VIEW SETMAGNIFICATION
          MAG = |MagCurrent|
      $!Endif
      $!If |DoZoom| == 1
        $!Varset |ActualMagFactor| = |MagFactor|
        $!If |SlowMo| != 0
          $!If |Loop| <= |SlowMo|
            $!Varset |ActualMagFactor| = (|ActualMagFactor|**0.5)
          $!Endif
          $!If |Loop| > (|NumSteps|-|SlowMo|)
            $!Varset |ActualMagFactor| = (|ActualMagFactor|**0.5)
          $!Endif
        $!Endif
        $!Varset |MagCurrent| = (|MagCurrent|*|ActualMagFactor|)
        $!VIEW SETMAGNIFICATION
          MAG = |MagCurrent|
      $!Endif
      ###
      ### Advance rotation
      $!If |DoRotate| == 1
        $!Varset |ActualRotation| = |RotationAngle|
        $!If |SlowMo| != 0
          $!Varset |ActualRotation| = (|RotationAngle| * |StepFactor|)
        $!Endif
        $!ROTATE3DVIEW |RotationAxis|
          ANGLE = |ActualRotation|
          ROTATEORIGINLOCATION = DEFINEDORIGIN
      $!Endif
      ###
      ### Redraw and save frame
      $!DRAWGRAPHICS TRUE
      $!If |Animation| >= 0
        $!Redraw
      $!Endif
      $!Varset |Counter| += 1
      $!If |Animation| == 1
        $!EXPORTNEXTFRAME 
      $!Endif
      $!If |Animation| > 1
        $!EXPORTSETUP EXPORTFNAME = "|ImagePath|/|AnimType||Counter%04d|.|Extension|"
        $!EXPORT
      $!Endif
      $!If |DoExtract| == 1
        $!RUNMACROFUNCTION "Extract Lines"
      $!Endif
    $!Endloop
  $!Endloop
  ###
  ### Finish up
  $!If |Animation| == -1
    $!Redraw
    $!Pause "Final view."
  $!Endif
  $!If |Animation| == 1
    $!EXPORTFINISH 
  $!Endif
  $!If |SaveRestoreStyle| == 1
    $!PromptForTextString |ReloadStyle|
      Instructions = "Animation is complete. Do you want to return to original view?\n0=NO  1=YES"
    $!If |ReloadStyle| == 1
      $!RUNMACROFUNCTION "Load Temporary Style"
    $!Endif
  $!Endif
$!Endmacrofunction


$!Macrofunction Name = "Extract Lines"
  ShowInMacroPanel = False
  ###
  $!Varset |TmpZone| = (|NUMZONES| + 1)
  $!CREATESTREAMZONES 
    CONCATENATE = NO
  $!WRITEDATASET  "Line|Counter%04d|.dat" 
    INCLUDETEXT = NO
    INCLUDEGEOM = NO
    INCLUDECUSTOMLABELS = NO
    ASSOCIATELAYOUTWITHDATAFILE = NO
    ZONELIST =  [|TmpZone|-|NUMZONES|]
    VARPOSITIONLIST =  [1-3]
    BINARY = NO
    USEPOINTFORMAT = YES
    PRECISION = 8
  $!DELETEZONES  [|TmpZone|-|NUMZONES|]
$!Endmacrofunction


$!Macrofunction Name = "Zone Point B traces"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!PROMPTFORTEXTSTRING |trackzone|
    INSTRUCTIONS = "Enter the zone to use for B traces."
  $!RUNMACROFUNCTION "Save Temporary Style"
  $!ACTIVEFIELDZONES = [|trackzone|]
  $!VARSET |nPoints| = (|MaxI|)
  $!RUNMACROFUNCTION "Load Temporary Style"
  $!LOOP |nPoints|
    $!VARSET |looppoint| = (|loop|)
    $!GetFieldValue |xloc|
      Zone = |trackzone|
      Var = 1
      Index = |looppoint|
    $!GetFieldValue |yloc|
      Zone = |trackzone|
      Var = 2
      Index = |looppoint|
    $!GetFieldValue |zloc|
      Zone = |trackzone|
      Var = 3
      Index = |looppoint|
    $!STREAMTRACE ADD
      STREAMTYPE = VOLUMELINE
      DIRECTION = BOTH
      STARTPOS
        {
        X = |xloc|
        Y = |yloc|
        Z = |zloc|
        }
  $!Endloop
  $!DRAWGRAPHICS TRUE
$!Endmacrofunction


$!Varset |tType| = 0
$!Varset |nLines| = 0
$!MacroFunction Name = "Auto U/B Traces"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!RUNMACROFUNCTION "Save Temporary Style"
  $!PAUSE "Make sure that the 3D zone is on and vector variables\nare set.  Assumes 3D zone is first zone!"
  ###  Read type of traces
  $!PromptForTextString |tType|
    Instructions = "Enter the trace type:\n1=streamlines  2=fieldlines-equator  3=fieldlines-pole"
  ###  Read lines to create
  $!PromptForTextString |nLines|
    Instructions = "Enter the number of evenly spaced traces"
  $!RUNMACROFUNCTION "UB Trace"
  $!RUNMACROFUNCTION "Load Temporary Style"
  $!DRAWGRAPHICS TRUE
$!Endmacrofunction


$!MacroFunction Name = "30 B Traces"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!Varset |tType| = 2
  $!Varset |nLines| = 30
  $!RUNMACROFUNCTION "UB Trace"
  $!DRAWGRAPHICS TRUE
$!Endmacrofunction


$!MacroFunction Name = "Save B For Surface"
  ShowInMacroPanel = True
  $!DRAWGRAPHICS FALSE
  $!RUNMACROFUNCTION "Save Temporary Style"
  $!Varset |tType| = 2
  $!Varset |nLines| = 72
  $!RUNMACROFUNCTION "UB Trace"
  $!Varset |tmpZone| = (|NUMZONES|-|nLines|+1)
  $!WRITEDATASET  "extracted_lines.dat" 
    INCLUDETEXT = NO
    INCLUDEGEOM = NO
    INCLUDECUSTOMLABELS = NO
    ASSOCIATELAYOUTWITHDATAFILE = NO
    ZONELIST =  [|tmpZone|-|NUMZONES|]
    VARPOSITIONLIST =  [1-3]
    BINARY = NO
    USEPOINTFORMAT = YES
    PRECISION = 9
  $!RUNMACROFUNCTION "Load Temporary Style"
  $!DRAWGRAPHICS TRUE
$!Endmacrofunction


$!MacroFunction Name = "UB Trace"
  ShowInMacroPanel = False
  ###  Check if variable tR is already there, if not create it and mark it
  ###  for deletion at script completion.  Then set the value.
  ###  If it already exists, stop script with error.
  $!Varset |vars0| = |NUMVARS|
  $!ALTERDATA 
    EQUATION = '{TR} = sqrt(V1**2 + V2**2 + V3**2)' 
  $!If |NUMVARS| == |vars0|
    $!PromptForTextString |error|
      Instructions = "You must cancel script, variable TR already exists."
  $!Endif
  $!Varset |varTR| = |NUMVARS|
  ###  Set some additional values before beginning loop
  $!Varset |rTol| = 0.001
  $!Varset |WriteEach| = 0
  $!Varset |maxstreamsteps| = 3000
  $!GLOBALSTREAM MAXSTEPS = |maxstreamsteps|
  $!GLOBALSTREAM CELLFRACTION = 0.25
  $!GLOBALSTREAM MINCELLFRACTION = 1E-05
  $!Varset |xloc| = 0.0
  $!Varset |yloc| = 0.0
  $!Varset |zloc| = 0.0
  $!ACTIVEFIELDZONES = [1]
  $!Varset |xmin3d| = |minV1|
  $!Varset |ymin3d| = |minV2|
  $!Varset |zmin3d| = |minV3|
  $!Varset |xmax3d| = |maxV1|
  $!Varset |ymax3d| = |maxV2|
  $!Varset |zmax3d| = |maxV3|
  ###  Start loop to find traces
  $!Loop |nLines|
    $!Varset |xlockeep| = 0.0
    $!Varset |ylockeep| = 0.0
    $!Varset |zlockeep| = 0.0
    $!Varset |Angle| = (((|Loop|-1.)/|nLines|)*(2.*|PI|))
    $!Varset |Radius| = 8.
    $!Varset |delRadius| = 8.
    $!Varset |dir| = 1
    $!Varset |sinA| = (sin(|Angle|))
    $!Varset |cosA| = (cos(|Angle|))
    $!Varset |MaxLoops| = 50
    $!Varset |nLoops| = 0
    $!Varset |done| = 0
    $!Varset |error| = 0
    $!If |tType| == 1
      $!BLANKING VALUE{INCLUDE = YES}
      $!BLANKING VALUE{CONSTRAINT 1 {INCLUDE = YES}}
      $!BLANKING VALUE{CONSTRAINT 1 {VARA = 1}}
      $!BLANKING VALUE{CONSTRAINT 1 {VALUECUTOFF = -1.0}}
      $!BLANKING VALUE{CONSTRAINT 1 {CONSTRAINTOP2MODE = USECONSTANT}}
      $!BLANKING VALUE{CONSTRAINT 1 {RELOP = LESSTHANOREQUAL}}
    $!Endif
    $!CREATERECTANGULARZONE 
      IMAX = 1
      JMAX = 1
      KMAX = 1
      X1 = 0
      Y1 = 0
      Z1 = 0
      X2 = 0
      Y2 = 0
      Z2 = 0
      XVAR = 1
      YVAR = 2
      ZVAR = 3
    $!Varset |TextZone| = |NUMZONES|
    $!ACTIVEFIELDZONES = [|TextZone|]
    $!ATTACHTEXT 
      XYPOS
        {
        X = 10
        Y = 50
        }
      TEXTSHAPE
        {
        HEIGHT = 40
        }
      ZONE = |TextZone|
      ATTACHTOZONE = YES
      ANCHOR = HEADLEFT
      TEXT = 'Starting Line |loop| of |nLines|'
    $!DRAWGRAPHICS TRUE
    $!Redraw
    $!DRAWGRAPHICS FALSE
    $!DELETEZONES  [|TextZone|]
    $!While |done| == 0
      $!Varset |nLoops| += 1
      $!If |nLoops| == |MaxLoops|
        $!Varset |done| = 1
      $!Endif
      $!If |tType| == 1
        $!Varset |yloc| = (|Radius|*|sinA|)
        $!Varset |zloc| = (|Radius|*|cosA|)
      $!Endif
      $!If |tType| == 2
        $!Varset |xloc| = (|Radius|*|sinA|)
        $!Varset |yloc| = (|Radius|*|cosA|)
      $!Endif
      $!If |tType| == 3
        ### NOTE: Radius is really angle in degrees
        $!Varset |zloc| = (cos(|Radius|*|d2r|))
        $!Varset |xloc| = (sin(|Radius|*|d2r|)*|sinA|)
        $!Varset |yloc| = (sin(|Radius|*|d2r|)*|cosA|)
      $!Endif
      ###  Trace and analyze
      $!ACTIVEFIELDZONES = [1]
      $!STREAMTRACE ADD
        STREAMTYPE = VOLUMELINE
        DIRECTION = BOTH
        STARTPOS
          {
          X =       (|xloc|)
          Y =       (|yloc|)
          Z =       (|zloc|)
          }
      $!CREATESTREAMZONES 
          CONCATENATE = NO
      $!STREAMTRACE DELETEALL
      $!Varset |tmpZone| = |NUMZONES|
      $!ACTIVEFIELDZONES = [|tmpZone|]
      $!GETFIELDVALUE |xstart|
         ZONE = |tmpZone|
         VAR = 1
         INDEX = 1
      $!GETFIELDVALUE |ystart|
         ZONE = |tmpZone|
         VAR = 2
         INDEX = 1
      $!GETFIELDVALUE |zstart|
         ZONE = |tmpZone|
         VAR = 3
         INDEX = 1
      $!GETFIELDVALUE |rstart|
         ZONE = |tmpZone|
         VAR = |varTR|
         INDEX = 1
      $!GETFIELDVALUE |xend|
         ZONE = |tmpZone|
         VAR = 1
         INDEX = |maxi|
      $!GETFIELDVALUE |yend|
         ZONE = |tmpZone|
         VAR = 2
         INDEX = |maxi|
      $!GETFIELDVALUE |zend|
         ZONE = |tmpZone|
         VAR = 3
         INDEX = |maxi|
      $!GETFIELDVALUE |rend|
         ZONE = |tmpZone|
         VAR = |varTR|
         INDEX = |maxi|
      $!DELETEZONES  [|tmpZone|]
      $!Varset |value| = 1
      $!If |tType| == 1
        ### stream: if start or end is within 1 of max X, then hits upstream boundary
        $!If (abs(|xstart|-|xmax3d|)) < 1.
          $!Varset |value| = -1
        $!Endif
        $!If (abs(|xend|-|xmax3d|)) < 1.
          $!Varset |value| = -1
        $!Endif
        $!If |value| < 0
          $!Varset |xlockeep| = |xloc|
          $!Varset |ylockeep| = |yloc|
          $!Varset |zlockeep| = |zloc|
        $!Endif
      $!Endif
      $!If |tType| == 2
        ###  field: if both ends are within a radius of 3, then closed fieldline
        $!If |rstart| > 3.
          $!Varset |value| = -1
        $!Endif
        $!If |rend| > 3.
          $!Varset |value| = -1
        $!Endif
        $!If |value| > 0
          $!Varset |xlockeep| = |xloc|
          $!Varset |ylockeep| = |yloc|
          $!Varset |zlockeep| = |zloc|
        $!Endif
      $!Endif
      $!If |tType| == 3
        ###  field: if both ends are within a radius of 3, then closed fieldline
        $!Varset |value| = -1
        $!If |rstart| > 3.
          $!Varset |value| = 1
        $!Endif
        $!If |rend| > 3.
          $!Varset |value| = 1
        $!Endif
        $!If |value| < 0
          $!Varset |xlockeep| = |xloc|
          $!Varset |ylockeep| = |yloc|
          $!Varset |zlockeep| = |zloc|
        $!Endif
      $!Endif
      $!If |value| > 0
        $!If |dir| == -1
          $!Varset |delRadius| = (|delRadius|/(-2.5))
          $!Varset |dir| = 1
        $!Endif
      $!Endif
      $!If |value| < 0
        $!If |dir| == 1
          $!Varset |delRadius| = (|delRadius|/(-2.5))
          $!Varset |dir| = -1
        $!Endif
      $!Endif
      $!Varset |Radius| = (|Radius| + |delRadius|)
      $!If (abs(|delRadius|)) < |rTol|
        $!Varset |done| = 1
      $!Endif
      $!If |Radius| < 1.
        $!Varset |done| = 1
        $!Varset |error| = 1
      $!Endif
    $!Endwhile
    $!If |error| == 0
      ### Add stream and extract to zone
      $!BLANKING VALUE{INCLUDE = NO}
      $!ACTIVEFIELDZONES = [1]
      $!STREAMTRACE ADD
        STREAMTYPE = VOLUMELINE
        DIRECTION = BOTH
        STARTPOS
          {
          X = |xlockeep|
          Y = |ylockeep|
          Z = |zlockeep|
          }
      $!CREATESTREAMZONES 
          CONCATENATE = NO
      $!STREAMTRACE DELETEALL
    $!Endif
    $!If |WriteEach| == 1
      $!WRITEDATASET  "extracted_trace_|loop%03d|.dat" 
        INCLUDETEXT = NO
        INCLUDEGEOM = NO
        INCLUDECUSTOMLABELS = NO
        ASSOCIATELAYOUTWITHDATAFILE = NO
        ZONELIST =  [|NUMZONES|]
        VARPOSITIONLIST =  [1-3]
        BINARY = NO
        USEPOINTFORMAT = YES
        PRECISION = 9
    $!Endif
  $!Endloop
  $!DELETEVARS  [|varTR|]
$!EndMacroFunction


########################################################### PAGE 4
########################################################### PAGE 4
########################################################### PAGE 4


$!MacroFunction Name = "-=- Amtec Macros -=-"
$!EndMacroFunction


#
# This macro function is used by the Cascade and Tile
# frames macros. It will prompt the user for PaperWidth
# and PaperHeight dimensions.  If the values are invalid,
# it will require the user to retype the values.
#
$!MACROFUNCTION
  NAME = "GetPaperDim"
  SHOWINMACROPANEL = FALSE  #We don't want this to show in the quick macro panel
  #
  # If you always use one paper size, instead of prompting
  # for the paper size, you can just set it explicitly.  Example:
  #
   $!VARSET |PAPERWIDTH| = 11
   $!VARSET |PAPERHEIGHT| = 8.5
  #
  
  #
  # Get the paper width
  #
 # $!VARSET |PAPERDIMOK| = 0
 # $!WHILE |PAPERDIMOK| == 0
 #   $!PROMPTFORTEXTSTRING |PAPERWIDTH|
 #     INSTRUCTIONS = "Enter the paper width."
 #   $!VARSET |PAPERDIMOK| = 1
 #   $!IF |PAPERWIDTH| < 1
 #     $!VARSET |PAPERDIMOK| = 0
 #     $!PAUSE "Paper width must be greater than or equal to 1."
 #   $!ENDIF
 # $!ENDWHILE
  
  #
  # Get the paper height
  #
 # $!VARSET |PAPERDIMOK| = 0
 # $!WHILE |PAPERDIMOK| == 0
 #   $!PROMPTFORTEXTSTRING |PAPERHEIGHT|
 #     INSTRUCTIONS = "Enter the paper height."
 #   $!VARSET |PAPERDIMOK| = 1
 #   $!IF |PAPERHEIGHT| < 1
 #     $!VARSET |PAPERDIMOK| = 0
 #     $!PAUSE "Paper height must be greater than or equal to 1."
 #   $!ENDIF
 # $!ENDWHILE  
$!ENDMACROFUNCTION


#
# Function Cascade Frames - Scott Fowler 11/2001
#
# This macro function will move and resize all the frames
# the layout in a cascading fashion. If there are too many
# frames to fit in one cascade, multiple layers of cascading
# frames will be produced. Frame order is maintained.
#
$!MACROFUNCTION
  NAME = "Cascade Frames"

  $!RUNMACROFUNCTION "GetPaperDim"

  # Change the delta to change the distance 
  # in which the frame are cascaded
  $!VARSET |DELTA|  = 0.2

  # These two values define the margin between the paper
  # and the edge of the outer frames.
  $!VARSET |STARTX| = 0.5
  $!VARSET |STARTY| = 0.5

  # This calculates a width and height such that all frames are
  # the same size.
  $!VARSET |FRAMEDIMOK| = 0
  $!VARSET |NUMLAYERS| = 1
  $!VARSET |FRAMESPERLAYER| = |NUMFRAMES|
  $!WHILE |FRAMEDIMOK| == 0
    $!VARSET |FRAMEWIDTH|  = ((|PAPERWIDTH|)  - ((|FRAMESPERLAYER|-1)*|DELTA|) - (|STARTX|*2))
    $!VARSET |FRAMEHEIGHT| = ((|PAPERHEIGHT|) - ((|FRAMESPERLAYER|-1)*|DELTA|) - (|STARTY|*2))

    $!VARSET |FRAMEDIMOK| = 1
    $!IF |FRAMEWIDTH| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEHEIGHT| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEDIMOK| == 0
      $!VARSET |NUMLAYERS| += 1
      $!VARSET |FRAMESPERLAYER| = (ceil(|NUMFRAMES|/|NUMLAYERS|))
    $!ENDIF
  $!ENDWHILE

  #
  # Now, reposition and resize each frame.
  #
  $!DRAWGRAPHICS NO
  $!VARSET |NUMFRAMESDRAWN| = 0
  $!LOOP |NUMLAYERS|
    $!VARSET |XPOS| = |STARTX|
    $!VARSET |YPOS| = |STARTY|
    $!LOOP |FRAMESPERLAYER|
      $!IF |NUMFRAMESDRAWN| < |NUMFRAMES|
        $!FRAMECONTROL POP
          FRAME = 1
        $!IF |LOOP| != 1
          $!VARSET |XPOS| += |DELTA|
          $!VARSET |YPOS| += |DELTA|
        $!ENDIF
        $!FRAMELAYOUT XYPOS
          {
            X = |XPOS|
            Y = |YPOS|
          }
        $!FRAMELAYOUT WIDTH  = |FRAMEWIDTH|
        $!FRAMELAYOUT HEIGHT = |FRAMEHEIGHT|
        $!VARSET |NUMFRAMESDRAWN| += 1
      $!ENDIF
    $!ENDLOOP
  $!ENDLOOP
  $!DRAWGRAPHICS YES
$!ENDMACROFUNCTION


#
# Function Tile Frames - Scott Fowler 11/2001
#                        Revised 8/2002
#
# This macro function will move and resize all the frames
# the layout in a tile fashion.
#
# It takes one parameter - "DOHORIZONTAL" or "DOVERTICAL"
# depending on which way the frames are to be tiled.
#
# If there are too many frames using the specified |MARGIN|,
# the |MARGIN| will be reduced until all frames fit on the
# paper.  |MARGIN| can be negative which will result in 
# overlapping frames. Frame order is maintained.
#
$!MACROFUNCTION
  NAME = "Tile Frames"
  SHOWINMACROPANEL = FALSE

  $!RUNMACROFUNCTION "GetPaperDim"
  
  # Change the margin to change the distance 
  # between frames and the margin to the edge of the paper.
  # This value may be modified automatically if there are
  # too many frames for the specified area.
  $!VARSET |MARGIN|  = 0.0

  # These two values define the margin from the edge of the paper
  # to the edge of the outer frames.
  $!VARSET |XPOS|   = 0.5
  $!VARSET |YPOS|   = 0.5

  $!VARSET |NUMFRAMESACROSS| = (sqrt(|NUMFRAMES|))
  $!VARSET |TMP| = (int(|NUMFRAMESACROSS|))
  # If NUMFRAMESACROSS is not an integer, "round" up
  $!IF |NUMFRAMESACROSS| > |TMP|
    $!VARSET |NUMFRAMESACROSS| = (|TMP|)
    $!IF |PAPERWIDTH| >= |PAPERHEIGHT|
      $!VARSET |NUMFRAMESACROSS| += 1
    $!ENDIF
  $!ENDIF

  $!VARSET |NUMFRAMESDOWN| = (|NUMFRAMES|/|NUMFRAMESACROSS|)
  $!VARSET |TMP| = (int(|NUMFRAMESDOWN|))
  # If NUMFRAMESDOWN is not an integer, "round" up
  $!IF |NUMFRAMESDOWN| > |TMP|
    $!VARSET |NUMFRAMESDOWN| = (|TMP|+1)
  $!ENDIF

  #
  # If we want to tile horizontally, swap the
  # |NUMFRAMESACROSS| and |NUMFRAMESDOWN| since
  # the above calculations are for vertically tile
  # frames
  #
  $!IF "|1|" == "DOHORIZONTAL"
    $!VARSET |TMP|             = |NUMFRAMESACROSS|
    $!VARSET |NUMFRAMESACROSS| = |NUMFRAMESDOWN|
    $!VARSET |NUMFRAMESDOWN|   = |TMP|
  $!ENDIF

  # This calculates a width and height such that all frames are
  # the same size, and fit with the specified margin.
  

  $!VARSET |FRAMEDIMOK| = 0
  $!WHILE |FRAMEDIMOK| == 0
    $!VARSET |FRAMEWIDTH|  = ((|PAPERWIDTH|-(|XPOS|*2)-((|NUMFRAMESACROSS|-1)*|MARGIN|))/(|NUMFRAMESACROSS|))
    $!VARSET |FRAMEHEIGHT| = ((|PAPERHEIGHT|-(|YPOS|*2)-((|NUMFRAMESDOWN|-1)*|MARGIN|))/(|NUMFRAMESDOWN|))

    $!VARSET |FRAMEDIMOK| = 1
    $!IF |FRAMEWIDTH| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEHEIGHT| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEDIMOK| == 0
      $!VARSET |MARGIN| -= 0.01
    $!ENDIF
  $!ENDWHILE
 
  $!DRAWGRAPHICS NO

  $!VARSET |FRAMESDRAWN| = 1
  $!VARSET |NUMDOWNDRAWN| = 1

  $!LOOP |NUMFRAMESDOWN|        
    $!LOOP |NUMFRAMESACROSS|
      # Make sure that we want to draw this frame
      $!IF |FRAMESDRAWN| <= |NUMFRAMES|        
        $!FRAMECONTROL POP
          FRAME = 1
        $!VARSET |FRAMEX| = ((|LOOP|-1)*(|FRAMEWIDTH|+|MARGIN|)+|XPOS|)
        $!VARSET |FRAMEY| = ((|NUMDOWNDRAWN|-1)*(|FRAMEHEIGHT|+|MARGIN|)+|YPOS|)
        $!FRAMELAYOUT XYPOS
          {
            X = |FRAMEX|
            Y = |FRAMEY|
          }
        $!FRAMELAYOUT WIDTH  = |FRAMEWIDTH|
        $!FRAMELAYOUT HEIGHT = |FRAMEHEIGHT|
        $!VARSET |FRAMESDRAWN| += 1
      $!ENDIF
    $!ENDLOOP    
    # Increment the number of down rows.
    $!VARSET |NUMDOWNDRAWN| += 1
  $!ENDLOOP

  $!DRAWGRAPHICS YES
$!ENDMACROFUNCTION


$!MACROFUNCTION
  NAME = "Tile Frames Vertically"
  $!RUNMACROFUNCTION "Tile Frames" ("DOVERTICAL")
$!ENDMACROFUNCTION


$!MACROFUNCTION
  NAME = "Tile Frames Horizontally"
  $!RUNMACROFUNCTION "Tile Frames" ("DOHORIZONTAL")
$!ENDMACROFUNCTION


###
###
###

