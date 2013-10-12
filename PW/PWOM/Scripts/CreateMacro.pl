#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help = ($h or $help);
my $Debug = ($d or $debug);
my $Verbose = ($v or $verbose);

use strict;
use FileHandle;
use Cwd;

my $dir = cwd;
&print_help if $Help;

my $plot;

print STDOUT "Enter Alt Slice For Output: ";
my $AltSlice = <STDIN>;
chop($AltSlice);
my @plot_file = sort glob "plot_snapshot????_iAlt${AltSlice}.dat";

my $nfile = @plot_file;

open OUT, ">Macro.mcr";

my $iPlot=1;
foreach $plot (@plot_file) {
    if ($iPlot eq 1){
	
	print OUT "#!MC 1000
\$!VarSet |MFBD| = '${dir}' 
\$!READDATASET  '\"|MFBD|/${plot}\" ' 
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  INITIALPLOTFIRSTZONEONLY = YES
  VARLOADMODE = BYNAME
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '\"X\" \"Y\" \"Z\" \"uO\" \"uHe\" \"uH\" \"ue\" \"lgnO\" \"lgnHe\" \"lgnH\" \"lgne\" \"TO\" \"THe\" \"TH\" \"Te\" \"MO\" \"MH\" \"MHe\" \"Me\" \"Ef\" \"Pe\"' 
\$!FIELDLAYERS SHOWMESH = NO
\$!FIELDLAYERS SHOWBOUNDARY = NO
\$!TRIANGULATE 
  SOURCEZONES =  [1]
  USEBOUNDARY = NO
  INCLUDEBOUNDARYPTS = NO
  TRIANGLEKEEPFACTOR = 0.25
\$!DELETEZONES  [1]
\$!GLOBALCONTOUR 1  VAR = 8
\$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
\$!FIELDLAYERS SHOWCONTOUR = YES";
	
    }
    $iPlot++;
    print OUT "
\$!READDATASET  '\"|MFBD|/${plot}\" ' 
   READDATAOPTION = APPEND
   RESETSTYLE = NO
   INCLUDETEXT = NO
   INCLUDEGEOM = NO
   INCLUDECUSTOMLABELS = NO
   INITIALPLOTFIRSTZONEONLY = YES
   VARLOADMODE = BYNAME
   INITIALPLOTTYPE = CARTESIAN2D
   VARNAMELIST = '\"X\" \"Y\" \"Z\" \"uO\" \"uHe\" \"uH\" \"ue\" \"lgnO\" \"lgnHe\" \"lgnH\" \"lgne\" \"TO\" \"THe\" \"TH\" \"Te\" \"MO\" \"MH\" \"MHe\" \"Me\" \"Ef\" \"Pe\"' 
\$!TRIANGULATE 
   SOURCEZONES =  [${iPlot}]
   USEBOUNDARY = NO
   INCLUDEBOUNDARYPTS = NO
   TRIANGLEKEEPFACTOR = 0.25
\$!DELETEZONES  [${iPlot}]";
}

#add grid
print OUT "
\$!VIEW FIT
\$!ROTATE2DDATA 
  ANGLE = 90
  X = 0
  Y = 0
  ZONELIST =  [1-${iPlot}]
\$!REDRAW 
\$!ATTACHGEOM 
  GEOMTYPE = CIRCLE
  LINETHICKNESS = 0.4
  ANCHORPOS
    {
    X = 0.0
    Y = 0.0
    }
  RAWDATA
.173648177667 
\$!ATTACHGEOM 
  GEOMTYPE = CIRCLE
  LINETHICKNESS = 0.4
  ANCHORPOS
    {
    X = 0.0
    Y = 0.0
    }
  RAWDATA
.342020143326 
\$!ATTACHGEOM 
  GEOMTYPE = CIRCLE
  LINETHICKNESS = 0.4
  ANCHORPOS
    {
    X = 0.0
    Y = 0.0
    }
  RAWDATA
.5
\$!ATTACHGEOM 
  GEOMTYPE = CIRCLE
  LINETHICKNESS = 0.4
  ANCHORPOS
    {
    X = 0.0
    Y = 0.0
    }
  RAWDATA
.642787609687
";

exit;
