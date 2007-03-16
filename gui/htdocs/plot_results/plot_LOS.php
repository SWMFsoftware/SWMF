<?php    // Set posted value
$plotextension = ".plt";
$macroextension = ".mcr";

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function form1() {   //
  global $runpath, $runname, $cmp, $plottype, $plotstyle, $plotfile, $plotfilelist;
  global $PlotApplication, $countfiles;

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

  // Start form
  echo "
<FORM METHOD=\"post\" TARGET=\"$plotfile\" ACTION=\"createplot.php\">
<INPUT TYPE=hidden name=plotapplication value=\"$PlotApplication\">
<INPUT TYPE=hidden name=cmp value=\"$cmp\">
<INPUT TYPE=hidden name=plottype value=\"$plottype\">
<INPUT TYPE=hidden name=plotfile value=\"$plotfile\">
<INPUT TYPE=hidden name=filedir value=\"$runpath/$cmp\">
  "; 
  // CONTOUR
  $defcontMM = " CHECKED"; $defcontCUSTOM = "";
  if($defcont == "CUSTOM") {$defcontMM = ""; $defcontCUSTOM = " CHECKED";}
  $defconttypeV = " CHECKED"; $defconttypeD = ""; $defconttypeR = "";
  if($defconttype == "D") {$defconttypeV = ""; $defconttypeD = " CHECKED";}
  if($defconttype == "R") {$defconttypeV = ""; $defconttypeR = " CHECKED";}
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
<br>
&nbsp&nbsp&nbsp&nbsp 
<INPUT TYPE=radio name=conttype value=\"V\"$defconttypeV>Value [P]<br>
&nbsp&nbsp&nbsp&nbsp 
<INPUT TYPE=radio name=conttype value=\"D\"$defconttypeD>Difference [P2-P] <br>
&nbsp&nbsp&nbsp&nbsp 
<INPUT TYPE=radio name=conttype value=\"R\"$defconttypeR>Ratio [P2/(P+1e-10)] <br>
&nbsp&nbsp&nbsp&nbsp 
If needed, P2 = 
  ";
  $choice = "";
  if (eregi("_c2_", $plotfile)) { $choice = "c2"; }
  if (eregi("_c3_", $plotfile)) { $choice = "c3"; }
  $plotfilelist2 = GetXYZPlotList("$runpath/$cmp", "$choice");
  PrintListMenu($plotfilelist2, $defplotfile2, "plotfile2", "name");
  echo "
</div>
<INPUT TYPE=hidden name=plotview value=\"$choice\">
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

  $parameter = array('variable', 'cont', 'contmin', 'contmax', 'conttype', 'plotfile2', 'isgrid', 'textlabel', 'plotview');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  if ($cont == "CUSTOM") { $contcustom = ""; } else { $contcustom = "#"; }
  $conttypeV = "#"; $conttypeD = "#"; $conttypeR = "#";
  if ($conttype == "V") { $conttypeV = ""; }
  if ($conttype == "D") { $conttypeD = ""; }
  if ($conttype == "R") { $conttypeR = ""; }
  if ($isgrid == "YES") { $grid = ""; } else { $grid = "#"; }
  if ($textlabel) { $text = ""; } else { $text = "#"; }
  $plotviewc2 = "#"; $plotviewc3 = "#";
  if ($plotview == 'c2') { $plotviewc2 = ""; }
  if ($plotview == 'c3') { $plotviewc3 = ""; }
  Exec("mkdir $batchdir/$tmpdir;
        rsync -av $imagedir/style.sty $imagedir/script.mcr $imagedir/tecplot.mcr $batchdir/$tmpdir/;
        cd $batchdir/$tmpdir;
        cat script.mcr | 
sed s/REPLACEISC2/'$plotviewc2'/ | 
sed s/REPLACEISC3/'$plotviewc3'/ | 
sed s/REPLACEVAR/'$variable'/ | 
sed s/REPLACEISGRID/'$grid'/ | 
sed s/REPLACEISTEXT/'$text'/ | 
sed s/REPLACETEXT/'$textlabel'/ | 
sed s/REPLACEISCONTCUSTOM/'$contcustom'/ | 
sed s/REPLACECONTMIN/'$contmin'/ | 
sed s/REPLACECONTMIN/'$contmin'/ | 
sed s/REPLACECONTMAX/'$contmax'/ | 
sed s/REPLACECONTTYPEV/'$conttypeV'/ | 
sed s/REPLACECONTTYPED/'$conttypeD'/ | 
sed s/REPLACECONTTYPER/'$conttypeR'/ > batch.mcr");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot2() {   // Save defaults file
  global $imagedir, $number;

  $parameter = array('variable', 'cont', 'contmin', 'contmax', 'conttype', 'plotfile2', 'isgrid', 'textlabel', 'plotview');

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

  $parameter = array('plotfile2');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  Exec("cd $batchdir/$tmpdir;
        $tecplot -p batch.mcr -b ../../../$cmp/$plotfile ../../../$cmp/$plotfile2;
        cp -f batch.log ../.;
        mv print.cps $file1");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot4() {   // Download for self-plot instructions
  global $filedir, $plotfile, $imagedir, $number;

  $parameter = array('plotfile2');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  echo "
<h3>Download files to recreate this figure manually.</h3>
<A HREF=\"$filedir/$plotfile\">$plotfile</A><br>
<A HREF=\"$filedir/$plotfile2\">$plotfile2</A><br>
<A HREF=\"$imagedir/batch-${number}.mcr\">batch-${number}.mcr</A><br>
<A HREF=\"$imagedir/style.sty\">style.sty</A><br>
<A HREF=\"$imagedir/tecplot.mcr\">tecplot.mcr</A><br>
<br>
Then run<br>
<CENTER>tecplot -p batch-${number}.mcr $plotfile $plotfile2</CENTER>
  ";
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
?>
