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
  //  PrintPltForm($runname, $plotfilelist, $plotfile, $plotstyle);
  PrintXYZPltForm($runname, $plotfile, $plotstyle);
  echo "</div>\n";

  // Return if no plotfile
  if (! $plotfile) { return; }

  $variables = GetPlotVariables("$runpath/$cmp/$plotfile");
  $uvar = GetVectorVar("U", $variables);
  $bvar = GetVectorVar("B", $variables);
  $jvar = GetVectorVar("J", $variables);

  // Determine slice direction
  $slice = "";
  if(eregi("x=", $plotfile)) { $slice = "X"; }
  if(eregi("y=", $plotfile)) { $slice = "Y"; }
  if(eregi("z=", $plotfile)) { $slice = "Z"; }

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
<INPUT TYPE=hidden name=slice value=\"$slice\">
  "; 
  // CONTOUR
  $defcontMM = " CHECKED"; $defcontCUSTOM = "";
  if($defcont == "CUSTOM") {$defcontMM = ""; $defcontCUSTOM = " CHECKED";}
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
  // GRID
  $defisgridNO = " CHECKED"; $defisgridYES = "";
  if($defisgrid == "YES") {$defisgridNO = ""; $defisgridYES = " CHECKED";}
  echo "
<br><b>Grid:</b><br>
<div class=\"indent\">
Plot grid?
<INPUT TYPE=radio name=isgrid value=\"NO\"$defisgridNO>No
<INPUT TYPE=radio name=isgrid value=\"YES\"$defisgridYES>Yes
</div>
  "; 
  // VIEW
  echo "
<br><b>View:</b><br>
<div class=\"indent\">
Center at:
X=<INPUT type=text size=3 name=viewx value=\"$defviewx\">
Y=<INPUT type=text size=3 name=viewy value=\"$defviewy\">
with view width <INPUT type=text size=4 name=vieww value=\"$defvieww\">
</div>
  "; 
  // VECTOR TRACES
  $defisblinesNO = " CHECKED"; $defisblinesYES = "";
  if($defisblines == "YES") {$defisblinesNO = ""; $defisblinesYES = " CHECKED";}
  $defstreamisBLACK = " CHECKED"; $defstreamisWHITE = "";
  if($defstreamcolor == "WHITE") {$defstreamisBLACK = ""; $defstreamisWHITE = " CHECKED";}
  echo "
<br><b>Vector Traces:</b><br>
<div class=\"indent\">
Plot fieldlines?
<INPUT TYPE=radio name=isblines value=\"NO\"$defisblinesNO>No
<INPUT TYPE=radio name=isblines value=\"YES\"$defisblinesYES>Yes
&nbsp&nbsp Line Color:
<INPUT TYPE=radio name=streamcolor value=\"BLACK\"$defstreamisBLACK>Black
<INPUT TYPE=radio name=streamcolor value=\"WHITE\"$defstreamisWHITE>White
</div>
  "; 
  // BODY
  $defiscircleNO = " CHECKED"; $defiscircleYES = "";
  if($defiscircle == "YES") {$defiscircleNO = ""; $defiscircleYES = " CHECKED";}
  echo "
<br><b>Body:</b><br>
<div class=\"indent\">
Plot circle at origin?
<INPUT TYPE=radio name=iscircle value=\"NO\"$defiscircleNO>No
<INPUT TYPE=radio name=iscircle value=\"YES\"$defiscircleYES>Yes
with radius <INPUT type=text size=4 name=circler value=\"$defcircler\">
</div>
  "; 
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

  $parameter = array('variable', 'uvar', 'bvar', 'jvar', 'cont', 'contmin', 'contmax', 'slice', 'isgrid', 'viewx', 'viewy', 'vieww', 'isblines', 'streamcolor', 'iscircle', 'circler', 'textlabel');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  if ($cont == "CUSTOM") { $contcustom = ""; } else { $contcustom = "#"; }
  if ($isblines == "YES") { $vector = ""; } else { $vector = "#"; }
  if ($slice == "X") { $vector1 = $bvar+1; $vector2 = $bvar+2; }
  if ($slice == "Y") { $vector1 = $bvar  ; $vector2 = $bvar+2; }
  if ($slice == "Z") { $vector1 = $bvar  ; $vector2 = $bvar+1; }
  if ($iscircle == "YES") { $circle = ""; } else { $circle = "#"; }
  if ($isgrid == "YES") { $grid = ""; } else { $grid = "#"; }
  if ($textlabel) { $text = ""; } else { $text = "#"; }
  $viewx1 = $viewx-$vieww/2;
  $viewx2 = $viewx+$vieww/2;
  $viewy1 = $viewy-$vieww/2;
  $viewy2 = $viewy+$vieww/2;
  Exec("mkdir $batchdir/$tmpdir;
        rsync -av $imagedir/style.sty $imagedir/script.mcr $imagedir/tecplot.mcr $batchdir/$tmpdir/;
        cd $batchdir/$tmpdir;
        cat script.mcr | 
sed s/REPLACEVAR/'$variable'/ | 
sed s/REPLACECIRCLER/'$circler'/ | 
sed s/REPLACEISCIRCLE/'$circle'/ | 
sed s/REPLACEISGRID/'$grid'/ | 
sed s/REPLACEISBVECTOR/'$vector'/ | 
sed s/REPLACEVECTOR1/'$vector1'/ | 
sed s/REPLACEVECTOR2/'$vector2'/ | 
sed s/REPLACESTREAMCOLOR/'$streamcolor'/ | 
sed s/REPLACEISTEXT/'$text'/ | 
sed s/REPLACETEXT/'$textlabel'/ | 
sed s/REPLACEVIEWX1/'$viewx1'/ | 
sed s/REPLACEVIEWX2/'$viewx2'/ | 
sed s/REPLACEVIEWY1/'$viewy1'/ | 
sed s/REPLACEVIEWY2/'$viewy2'/ | 
sed s/REPLACEISCONTCUSTOM/'$contcustom'/ | 
sed s/REPLACECONTMIN/'$contmin'/ | 
sed s/REPLACECONTMIN/'$contmin'/ | 
sed s/REPLACECONTMAX/'$contmax'/ | 
sed s/REPLACESLICE/'$slice'/ > batch.mcr");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot2() {   // Save defaults file
  global $imagedir, $number;

  $parameter = array('variable', 'uvar', 'bvar', 'jvar', 'cont', 'contmin', 'contmax', 'slice', 'isgrid', 'viewx', 'viewy', 'vieww', 'isblines', 'streamcolor', 'iscircle', 'circler', 'textlabel');

  $deffile = "$imagedir/defaults-${number}.php";
  if(!$fh = fopen($deffile, 'w')) {
    echo "<h1>The file $deffile cannot be opened</h1>";
  } else {
    fwrite($fh, "<?php\n");
    foreach($parameter as $name) {
      $value = $_POST[$name];
      fwrite($fh, "\$def$name = \"$value\";\n");
    }
    $value = $_POST['slice'];
    fwrite($fh, "\$choice = \"$value\";\n");
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
        echo '<A HREF=\"${imagedir}/tecplot.mcr\">tecplot.mcr</A><br>' >> page.html;
        echo '<br>' >> page.html;
        echo 'Then run<br>' >> page.html;
        echo '<CENTER>tecplot -p batch-${number}.mcr ${plotfile}</CENTER>' >> page.html;
        echo '' >> page.html");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
?>
