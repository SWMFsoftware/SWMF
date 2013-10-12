<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function form1() {   //
  global $runpath, $runname, $cmp, $plottype, $plotstyle, $plotfile, $plotfilelist;
  global $countfiles;

  include("$runpath/images/${cmp}_$plottype/defaults-$plotstyle.php");

  echo "
<b>Plotfile:</b> T=Hour:Min:Sec N=Iterations ($countfiles files found)<br>
<div class=\"indent\">
  ";
  PrintPltForm($runname, $plotfilelist, $plotfile, $plotstyle);
  echo "</div>\n";

  // Return if no plotfile
  if (! $plotfile) { return; }

  $variables = GetPlotVariables("$runpath/$cmp/$plotfile");
  $uvar = GetVectorVar("U", $variables);
  $bvar = GetVectorVar("B", $variables);
  $jvar = GetVectorVar("J", $variables);

  // Set defaults
  $defcontMM = " CHECKED"; $defcontCUSTOM = "";
  if($defcont == "CUSTOM") {$defcontMM = ""; $defcontCUSTOM = " CHECKED";}
  $defsliceNO = ""; $defsliceX = ""; $defsliceY = " CHECKED"; $defsliceZ = "";
  if($defslice == "NO") {$defsliceNO = " CHECKED"; $defsliceX = ""; $defsliceY = ""; $defsliceZ = "";}
  if($defslice == "X")  {$defsliceNO = ""; $defsliceX = " CHECKED"; $defsliceY = ""; $defsliceZ = "";}
  if($defslice == "Z")  {$defsliceNO = ""; $defsliceX = ""; $defsliceY = ""; $defsliceZ = " CHECKED";}
  $defslice2NO = ""; $defslice2X = ""; $defslice2Y = " CHECKED"; $defslice2Z = "";
  if($defslice2 == "NO") {$defslice2NO = " CHECKED"; $defslice2X = ""; $defslice2Y = ""; $defslice2Z = "";}
  if($defslice2 == "X")  {$defslice2NO = ""; $defslice2X = " CHECKED"; $defslice2Y = ""; $defslice2Z = "";}
  if($defslice2 == "Z")  {$defslice2NO = ""; $defslice2X = ""; $defslice2Y = ""; $defslice2Z = " CHECKED";}
  $defisgridNO = " CHECKED"; $defisgridYES = "";
  if($defisgrid == "YES") {$defisgridNO = ""; $defisgridYES = " CHECKED";}
  $defisisoNO = " CHECKED"; $defisisoYES = "";
  if($defisiso == "YES") {$defisisoNO = ""; $defisisoYES = " CHECKED";}
  $defisbvectoreqNO = " CHECKED"; $defisbvectoreqYES = "";
  if($defisbvectoreq == "YES") {$defisbvectoreqNO = ""; $defisbvectoreqYES = " CHECKED";}
  $defissphereNO = " CHECKED"; $defissphereYES = "";
  if($defissphere == "YES") {$defissphereNO = ""; $defissphereYES = " CHECKED";}

  // Start form
  echo "
<FORM METHOD=\"post\" TARGET=\"$plotfile\" ACTION=\"createplot.php\">
<INPUT TYPE=hidden name=cmp value=\"$cmp\">
<INPUT TYPE=hidden name=plottype value=\"$plottype\">
<INPUT TYPE=hidden name=plotfile value=\"$plotfile\">
<INPUT TYPE=hidden name=filedir value=\"$runpath/$cmp\">
<INPUT TYPE=hidden name=uvar value=\"$uvar\">
<INPUT TYPE=hidden name=bvar value=\"$bvar\">
<INPUT TYPE=hidden name=jvar value=\"$jvar\">
  "; 
  // CONTOUR
  echo "
<b>Contour</b><br>
<div class=\"indent\">
Variable: 
  ";
  PrintListMenu($variables, $defvariable, "variable", "number");
  echo "
Range:
<INPUT TYPE=radio name=cont value=\"MM\"$defcontMM>Min/Max
<INPUT TYPE=radio name=cont value=\"CUSTOM\"$defcontCUSTOM>Custom
<INPUT type=text size=5 name=contmin value=\"$defcontmin\">
<INPUT type=text size=5 name=contmax value=\"$defcontmax\">
</div>
  "; 
  // SLICE & GRID
  if($cmp == "IM") {
    echo "
<br><b>Grid:</b><br>
<div class=\"indent\">
Plot grid?
<INPUT TYPE=radio name=isgrid value=\"NO\"$defisgridNO>No
<INPUT TYPE=radio name=isgrid value=\"YES\"$defisgridYES>Yes
</div>
    "; 
  } else {
    echo "
<br><b>Slice:</b><br>
<div class=\"indent\">
Slice 1: 
<INPUT TYPE=radio name=slice value=\"NO\"$defsliceNO>No Slice
<INPUT TYPE=radio name=slice value=\"X\"$defsliceX>X=
<INPUT TYPE=radio name=slice value=\"Y\"$defsliceY>Y=
<INPUT TYPE=radio name=slice value=\"Z\"$defsliceZ>Z=
<INPUT type=text size=4 name=sliceval value=\"$defsliceval\">
 <br>
Slice 2: 
<INPUT TYPE=radio name=slice2 value=\"NO\"$defslice2NO>No Slice
<INPUT TYPE=radio name=slice2 value=\"X\"$defslice2X>X=
<INPUT TYPE=radio name=slice2 value=\"Y\"$defslice2Y>Y=
<INPUT TYPE=radio name=slice2 value=\"Z\"$defslice2Z>Z=
<INPUT type=text size=4 name=slice2val value=\"$defslice2val\">
 <br>
Plot grid on slices?
<INPUT TYPE=radio name=isgrid value=\"NO\"$defisgridNO>No
<INPUT TYPE=radio name=isgrid value=\"YES\"$defisgridYES>Yes
</div>
    "; 
  }
  // ISOSURFACE
  if($cmp != "IM") {
    echo "
<br><b>Isosurface:</b><br>
<div class=\"indent\">
<INPUT TYPE=radio name=isiso value=\"NO\"$defisisoNO>No
<INPUT TYPE=radio name=isiso value=\"YES\"$defisisoYES>Yes &nbsp &nbsp
Variable=
    ";
    PrintIsoVarsMenu($variables, $defisovar);
    echo "
Value=<INPUT type=text size=4 name=isoval value=\"$defisoval\">
</div>
    "; 
  }
  // VIEW
  if($cmp != "IM") {
    echo "
<br><b>View:</b><br>
<div class=\"indent\">
Center at:
X=<INPUT type=text size=3 name=viewx value=\"$defviewx\">
Y=<INPUT type=text size=3 name=viewy value=\"$defviewy\">
Z=<INPUT type=text size=3 name=viewz value=\"$defviewz\">
with view width <INPUT type=text size=4 name=vieww value=\"$defvieww\">
 <br>
Perspective angles: 
Phi=<INPUT type=text size=4 name=viewphi value=\"$defviewphi\">
Theta=<INPUT type=text size=4 name=viewtheta value=\"$defviewtheta\">
(<a href=\"viewperspective.php\" TARGET=\"_angles\">Help me with view angles.</a>)
</div>
    "; 
  }
  // VECTOR TRACES
  if($cmp != "IM") {
    echo "
<br><b>Vector Traces:</b><br>
<div class=\"indent\">
Plot last closed fieldlines?
<INPUT TYPE=radio name=isbvectoreq value=\"NO\"$defisbvectoreqNO>No
<INPUT TYPE=radio name=isbvectoreq value=\"YES\"$defisbvectoreqYES>Yes  
(5-10 minutes render time)
</div>
    "; 
  }
  // BODY
  if($cmp != "IM") {
    echo "
<br><b>Body:</b><br>
<div class=\"indent\">
Plot sphere at origin?
<INPUT TYPE=radio name=issphere value=\"NO\"$defissphereNO>No
<INPUT TYPE=radio name=issphere value=\"YES\"$defissphereYES>Yes
with radius <INPUT type=text size=4 name=spherer value=\"$defspherer\">
</div>
    "; 
  }
  // TEXT LABEL
  echo "
<br><b>Text Label:</b><br>
<div class=\"indent\">
Label:
<INPUT type=text size=40 name=textlabel value=\"$deftextlabel\">
</div>
 <br>
<INPUT TYPE=submit VALUE=\"Update Plot\">   (wait ~1 minute unless fieldlines plotted)</FORM>
  ";
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot0() {   // Check for errors
  $parameter = array('variable');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  if ($variable) { return 0; }
  return 1;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot1() {   // Create scripts
  global $batchdir, $tmpdir, $imagedir;

  $parameter = array('variable', 'uvar', 'bvar', 'jvar', 'cont', 'contmin', 'contmax', 'slice', 'sliceval', 'slice2', 'slice2val', 'isgrid', 'isiso', 'isovar', 'isoval', 'viewx', 'viewy', 'viewz', 'vieww', 'viewphi', 'viewtheta', 'isbvectoreq', 'issphere', 'spherer', 'textlabel');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  if ($cont == "CUSTOM") { $contcustom = ""; } else { $contcustom = "#"; }
  if ($isbvectoreq == "YES") { $vector = ""; } else { $vector = "#"; }
  if ($isbvectoreq == "YES") { $vector1 = $bvar; $vector2 = $bvar+1; $vector3 = $bvar+2; }
    else                     { $vector1 = "";    $vector2 = "";      $vector3 = "";      }
  if ($issphere == "YES") { $sphere = ""; } else { $sphere = "#"; }
  if ($isgrid == "YES") { $grid = ""; } else { $grid = "#"; }
  if ($isiso == "YES") { $iso = "";  $notiso = "#"; }
    else               { $iso = "#"; $notiso = "";  }
  if ($slice == "X")   { $slicex  = "X=1"; } else { $slicex  = "X=0"; }
  if ($slice == "Y")   { $slicey  = "Y=1"; } else { $slicey  = "Y=0"; }
  if ($slice == "Z")   { $slicez  = "Z=1"; } else { $slicez  = "Z=0"; }
  if ($slice == "NO")  { $isnotslice1 = "";  $isslice1 = "#"; }
    else               { $isnotslice1 = "#"; $isslice1 = "";  }
  if ($slice2 == "X")  { $slice2x = "X=1"; } else { $slice2x = "X=0"; }
  if ($slice2 == "Y")  { $slice2y = "Y=1"; } else { $slice2y = "Y=0"; }
  if ($slice2 == "Z")  { $slice2z = "Z=1"; } else { $slice2z = "Z=0"; }
  if ($slice2 == "NO") { $slice2z = "Z=1"; }  // Just make Z slice
  if ($slice2 == "NO") { $isnotslice2 = "";  $isslice2 = "#"; }
    else               { $isnotslice2 = "#"; $isslice2 = "";  }
  if ($textlabel) { $text = ""; } else { $text = "#"; }
  Exec("mkdir $batchdir/$tmpdir;
        rsync -av $imagedir/plts/* $batchdir/$tmpdir/;
        rsync -av $imagedir/style.sty $imagedir/script.mcr $imagedir/tecplot.mcr $batchdir/$tmpdir/;
        cd $batchdir/$tmpdir;
        cat script.mcr | 
sed s/REPLACEVAR/'$variable'/ | 
sed s/REPLACESPHERER/'$spherer'/ | 
sed s/REPLACEISSPHERE/'$sphere'/ | 
sed s/REPLACEISGRID/'$grid'/ | 
sed s/REPLACEISISO/'$iso'/ | 
sed s/REPLACEISNOTISO/'$notiso'/ | 
sed s/REPLACEISOVAR/'$isovar'/ | 
sed s/REPLACEISOVAL/'$isoval'/ | 
sed s/REPLACEISBVECTOR/'$vector'/ | 
sed s/REPLACEVECTOR1/'$vector1'/ | 
sed s/REPLACEVECTOR2/'$vector2'/ | 
sed s/REPLACEVECTOR3/'$vector3'/ | 
sed s/REPLACEISTEXT/'$text'/ | 
sed s/REPLACETEXT/'$textlabel'/ | 
sed s/REPLACEVIEWPHI/'$viewphi'/ | 
sed s/REPLACEVIEWTHETA/'$viewtheta'/ | 
sed s/REPLACEVIEWTOX/'$viewx'/ | 
sed s/REPLACEVIEWTOY/'$viewy'/ | 
sed s/REPLACEVIEWTOZ/'$viewz'/ | 
sed s/REPLACEVIEWTOW/'$vieww'/ | 
sed s/REPLACEISCONTCUSTOM/'$contcustom'/ | 
sed s/REPLACECONTMIN/'$contmin'/ | 
sed s/REPLACECONTMIN/'$contmin'/ | 
sed s/REPLACECONTMAX/'$contmax'/ | 
sed s/REPLACEISSLICE1/'$isslice1'/ | 
sed s/REPLACEISNOTSLICE1/'$isnotslice1'/ | 
sed s/REPLACEISSLICE2/'$isslice2'/ | 
sed s/REPLACEISNOTSLICE2/'$isnotslice2'/ | 
sed s/REPLACESLICENORMALX/'$slicex'/ | 
sed s/REPLACESLICENORMALY/'$slicey'/ | 
sed s/REPLACESLICENORMALZ/'$slicez'/ | 
sed s/REPLACESLICEVALUE/'$slice=$sliceval'/ |
sed s/REPLACESLICE2NORMALX/'$slice2x'/ | 
sed s/REPLACESLICE2NORMALY/'$slice2y'/ | 
sed s/REPLACESLICE2NORMALZ/'$slice2z'/ | 
sed s/REPLACESLICE2VALUE/'$slice2=$slice2val'/ > batch.mcr");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot2() {   // Save defaults file
  global $imagedir, $number;

  $parameter = array('variable', 'uvar', 'bvar', 'jvar', 'cont', 'contmin', 'contmax', 'slice', 'sliceval', 'slice2', 'slice2val', 'isgrid', 'isiso', 'isovar', 'isoval', 'viewx', 'viewy', 'viewz', 'vieww', 'viewphi', 'viewtheta', 'isbvectoreq', 'issphere', 'spherer', 'textlabel');

  $deffile = "$imagedir/defaults-${number}.php";
  if(!$fh = fopen($deffile, 'w')) {
    echo "<h1>The file $deffile cannot be opened</h1>";
  } else {
    fwrite($fh, "<?php\n");
    foreach($parameter as $name) {
      $value = $_POST[$name];
      fwrite($fh, "\$def$name = \"$value\";\n");
    }
    fwrite($fh, " ?>\n");
    fclose($fh);
  }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot3() {   // Run Script
  global $batchdir, $tmpdir, $cmp, $plotfile, $file1;
  include("paths.php");

  Exec("cd $batchdir/$tmpdir;
        echo '#!/bin/sh' > batchscript.sh;
        echo '' >> batchscript.sh;
        echo '${tecplot} -p batch.mcr -b ../../../${cmp}/${plotfile}' >> batchscript.sh;
        echo 'cp -f batch.log ../.' >> batchscript.sh;
        echo 'mv print.cps ${file1}' >> batchscript.sh;
        echo '' >> batchscript.sh;
        chmod 755 batchscript.sh");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot4() {   // Download for self-plot instructions
  global $batchdir, $tmpdir, $filedir, $plotfile, $imagedir, $number;

  Exec("cd $batchdir/$tmpdir;
        echo '<h3>Download files to recreate this figure manually.</h3>' >> page.html;
        echo '<A HREF=\"${filedir}/${plotfile}\">${plotfile}</A><br>' >> page.html;
        echo '<A HREF=\"${imagedir}/batch-${number}.mcr\">batch-${number}.mcr</A><br>' >> page.html;
        echo '<A HREF=\"${imagedir}/style.sty\">style.sty</A><br>' >> page.html;
        echo '<A HREF=\"${imagedir}/tecplot.mcr\">tecplot.mcr</A><br>' >> page.html");
  $dir2 = opendir("$imagedir/plts");
  while( $file = readdir( $dir2 ) ) {
    if (eregi(".plt", $file)) {
      Exec("cd $batchdir/$tmpdir;
            echo '<A HREF=\"${imagedir}/plts/${file}\">${file}</A><br>' >> page.html");
    }
  }
  Exec("cd $batchdir/$tmpdir;
        echo '<br>' >> page.html;
        echo 'Then run<br>' >> page.html;
        echo '<CENTER>tecplot -p batch-${number}.mcr ${plotfile}</CENTER>' >> page.html;
        echo '' >> page.html");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
?>
