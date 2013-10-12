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
<b>Logfile:</b> ($countfiles files found)<br>
<div class=\"indent\">
  ";
  PrintPltForm($runname, $plotfilelist, $plotfile, $plotstyle);
  echo "</div>\n";

  // Return if no plotfile
  if (! $plotfile) { return; }

  $variables = GetPlotVariables("$runpath/$cmp/$plotfile");

  // Start form
  echo "
<FORM METHOD=\"post\" TARGET=\"$plotfile\" ACTION=\"createplot.php\">
<INPUT TYPE=hidden name=cmp value=\"$cmp\">
<INPUT TYPE=hidden name=plottype value=\"$plottype\">
<INPUT TYPE=hidden name=plotfile value=\"$plotfile\">
<INPUT TYPE=hidden name=filedir value=\"$runpath/$cmp\">
  "; 
  // X-AXIS
  $defxaxisMM = " CHECKED"; $defxaxisCUSTOM = "";
  if($defxaxis == "CUSTOM") {$defxaxisMM = ""; $defxaxisCUSTOM = " CHECKED";}
  echo "
<b>X Axis</b><br>
<div class=\"indent\">
Variable: 
  ";
  PrintListMenu($variables, $defxvar, "xvar", "number");
  echo "
Range:
<INPUT TYPE=radio name=xaxis value=\"MM\"$defxaxisMM>Min/Max
<INPUT TYPE=radio name=xaxis value=\"CUSTOM\"$defxaxisCUSTOM>Custom
<INPUT type=text size=5 name=xaxismin value=\"$defxaxismin\">
<INPUT type=text size=5 name=xaxismax value=\"$defxaxismax\">
</div>
  "; 
  // Y-AXIS
  $defyaxisMM = " CHECKED"; $defyaxisCUSTOM = "";
  if($defyaxis == "CUSTOM") {$defyaxisMM = ""; $defyaxisCUSTOM = " CHECKED";}
  echo "
<b>Y Axis</b><br>
<div class=\"indent\">
Variable: 
  ";
  PrintListMenu($variables, $defyvar, "yvar", "number");
  echo "
Range:
<INPUT TYPE=radio name=yaxis value=\"MM\"$defyaxisMM>Min/Max
<INPUT TYPE=radio name=yaxis value=\"CUSTOM\"$defyaxisCUSTOM>Custom
<INPUT type=text size=5 name=yaxismin value=\"$defyaxismin\">
<INPUT type=text size=5 name=yaxismax value=\"$defyaxismax\">
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
  return 0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot1() {   // Create scripts
  global $batchdir, $tmpdir, $imagedir;

  $parameter = array('yvar', 'yaxis', 'yaxismin', 'yaxismax', 'xvar', 'xaxis', 'xaxismin', 'xaxismax', 'textlabel');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  if ($xaxis == "CUSTOM") { $xaxiscustom = ""; } else { $xaxiscustom = "#"; }
  if ($yaxis == "CUSTOM") { $yaxiscustom = ""; } else { $yaxiscustom = "#"; }
  if ($textlabel) { $text = ""; } else { $text = "#"; }
  Exec("mkdir $batchdir/$tmpdir;
        rsync -av $imagedir/style.sty $imagedir/script.mcr $imagedir/tecplot.mcr $batchdir/$tmpdir/;
        cd $batchdir/$tmpdir;
        cat script.mcr | 
sed s/REPLACEXVAR/'$xvar'/ | 
sed s/REPLACEISXAXISCUSTOM/'$xaxiscustom'/ | 
sed s/REPLACEXAXISMIN/'$xaxismin'/ | 
sed s/REPLACEXAXISMIN/'$xaxismin'/ | 
sed s/REPLACEXAXISMAX/'$xaxismax'/ | 
sed s/REPLACEYVAR/'$yvar'/ | 
sed s/REPLACEISYAXISCUSTOM/'$yaxiscustom'/ | 
sed s/REPLACEYAXISMIN/'$yaxismin'/ | 
sed s/REPLACEYAXISMIN/'$yaxismin'/ | 
sed s/REPLACEYAXISMAX/'$yaxismax'/ | 
sed s/REPLACEISTEXT/'$text'/ | 
sed s/REPLACETEXT/'$textlabel'/ > batch.mcr");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot2() {   // Save defaults file
  global $imagedir, $number;

  $parameter = array('yvar', 'yaxis', 'yaxismin', 'yaxismax', 'xvar', 'xaxis', 'xaxismin', 'xaxismax', 'textlabel');

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
        echo '<A HREF=\"${imagedir}/tecplot.mcr\">tecplot.mcr</A><br>' >> page.html;
        echo '<br>' >> page.html;
        echo 'Then run<br>' >> page.html;
        echo '<CENTER>tecplot -p batch-${number}.mcr ${plotfile}</CENTER>' >> page.html;
        echo '' >> page.html");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
?>
