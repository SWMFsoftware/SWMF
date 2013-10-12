#!MC 1000
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

$!Varset |PI| = (2.*asin(1.))
$!VarSet |ContMin| =  -1.0
$!VarSet |ContMax| =   1.0
$!VarSet |ContNum| = 201
$!VarSet |NN| = 0
$!VarSet |DIR| = 'NONE'
$!VarSet |VAR| =   'none'


$!Macrofunction Name = "-=- Custom Menu -=-"
$!Endmacrofunction

$!MacroFunction Name = "Set Contours"
  ShowInMacroPanel = True
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
$!EndMacroFunction

$!MacroFunction Name = "Set rho"
  $!VarSet |VAR| =   'rho'
  $!GLOBALCONTOUR 1  VAR = 4
  $!VarSet |ContMin| =   0.0
  $!VarSet |ContMax| =  28.0
  $!VarSet |ContNum| = 201
  $!RUNMACROFUNCTION  "Set Contours" 
$!EndMacroFunction

$!MacroFunction Name = "Set ux"
  $!VarSet |VAR| =   'ux'
  $!GLOBALCONTOUR 1  VAR = 5
  $!VarSet |ContMin| = -400.0
  $!VarSet |ContMax| =  400.0
  $!VarSet |ContNum| = 201
  $!RUNMACROFUNCTION  "Set Contours" 
$!EndMacroFunction

$!MacroFunction Name = "Set p"
  $!VarSet |VAR| =   'p'
  $!GLOBALCONTOUR 1  VAR = 11
  $!VarSet |ContMin| =   0.0
  $!VarSet |ContMax| =   0.2
  $!VarSet |ContNum| = 201
  $!RUNMACROFUNCTION  "Set Contours" 
$!EndMacroFunction

$!MacroFunction Name = "Set jy"
  $!VarSet |VAR| =   'jy'
  $!GLOBALCONTOUR 1  VAR = 13
  $!VarSet |ContMin| =  -0.0005
  $!VarSet |ContMax| =   0.0005
  $!VarSet |ContNum| = 201
  $!RUNMACROFUNCTION  "Set Contours" 
$!EndMacroFunction

$!MacroFunction Name = "Export PS"
  $!EXPORTSETUP EXPORTFORMAT = PS
  $!PRINTSETUP PALETTE = COLOR
  $!EXPORTSETUP IMAGEWIDTH = 800
  $!EXPORTSETUP EXPORTFNAME = 'Slice|DIR||NN%04d||VAR|.ps'
  $!EXPORT 
    EXPORTREGION = CURRENTFRAME
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

$!MacroFunction Name = "Black Grids/Boundaries"
  ShowInMacroPanel = True
  $!FIELD [1-|Numzones|]  MESH{COLOR = BLACK}
  $!FIELD [1-|Numzones|]  BOUNDARY{COLOR = BLACK}
  $!REDRAW
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

