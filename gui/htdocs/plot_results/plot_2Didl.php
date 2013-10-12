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
<b>Plotfile:</b>  T=Hour:Min:Sec N=Iterations ($countfiles files found)<br>
<div class=\"indent\">
  ";
  //  PrintPltForm($runname, $plotfilelist, $plotfile, $plotstyle);
  PrintXYZPltForm($runname, $plotfile, $plotstyle);
  echo "</div>\n";

  // Return if no plotfile
  if (! $plotfile) { return; }

  $variables = GetPlotVariables("$runpath/$cmp/$plotfile");

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
<INPUT TYPE=hidden name=slice value=\"$slice\">
  "; 
  // CONTOUR
  echo "
<b>Contour</b><br>
<div class=\"indent\">
Variable: 
  ";
  PrintListMenu($variables, $defvariable1, "variable1", "name");
  $plotoptions = array('contbar', 'contfill', 'contour');
  PrintListMenu($plotoptions, $defconttype1, "conttype1", "name");
  $checked1 = " CHECKED"; $checked2 = "";
  if($defcontrange1 == "CUSTOM") { $checked1 = ""; $checked2 = " CHECKED"; }
  echo "
Range:
<INPUT TYPE=radio name=contrange1 value=\"AUTO\"$checked1>Auto
<INPUT TYPE=radio name=contrange1 value=\"CUSTOM\"$checked2>Custom
<INPUT type=text size=5 name=contmin1 value=\"$defcontmin1\">
<INPUT type=text size=5 name=contmax1 value=\"$defcontmax1\">
</div>
  ";
  // GRID
  $checked1 = ""; $checked2 = " CHECKED";
  if($defgrid == "grid") { $checked1 = " CHECKED"; $checked2 = ""; }
  echo "
<br><b>Grid</b><br>
<div class=\"indent\">
Plot the grid?: 
<INPUT TYPE=radio name=grid value=\"grid\"$checked1>Yes
<INPUT TYPE=radio name=grid value=\"\"$checked2>No
</div>
  "; 
  // VECTORS
  echo "
<br><b>Vectors:</b><br>
<div class=\"indent\">
<table>
<tr><td>Plot velocity vectors?</td><td>
  ";
  $plotoptions = array('NO', ' vector', ' stream', ' stream2');
  PrintListMenu($plotoptions, $defuvector, "uvector", "name");
  echo "
</td></tr>
<tr><td>Plot field vectors?</td><td>
  ";
  $plotoptions = array('NO', ' vector', ' stream', ' stream2');
  PrintListMenu($plotoptions, $defbvector, "bvector", "name");
  echo "
</td></tr>
<tr><td>Plot current vectors?</td><td>
  ";
  $plotoptions = array('NO', ' vector', ' stream', ' stream2');
  PrintListMenu($plotoptions, $defjvector, "jvector", "name");
  echo "
</td></tr>
</table>
</div>
  "; 
  // VIEW
  $checked1 = " CHECKED"; $checked2 = "";
  if($defview == "Custom") { $checked1 = ""; $checked2 = " CHECKED"; }
  echo "
<br><b>View:</b><br>
<div class=\"indent\">
<INPUT TYPE=radio name=view value=\"\"$checked1>All
<INPUT TYPE=radio name=view value=\"Custom\"$checked2>Custom:
xmin, xmax=
<INPUT type=text size=3 name=viewx1 value=\"$defviewx1\">
<INPUT type=text size=3 name=viewx2 value=\"$defviewx2\">
ymin, ymax=
<INPUT type=text size=3 name=viewy1 value=\"$defviewy1\">
<INPUT type=text size=3 name=viewy2 value=\"$defviewy2\">
</div>
  "; 

  echo "
 <br>
<INPUT TYPE=submit VALUE=\"Update Plot\"></FORM>
  ";
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot0() {   // Check for errors
  $parameter = array('variable1');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  if ($variable1) { return 0; }
  return 1;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot1() {   // Create scripts
  global $imagedir, $filedir, $plotfile, $batchdir, $tmpdir, $macroextension;

  $parameter = array('variable1', 'conttype1', 'contrange1', 'contmin1', 'contmax1', 'slice', 'grid', 'uvector', 'bvector', 'jvector', 'view', 'viewx1', 'viewx2', 'viewy1', 'viewy2');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  $runidlexists = "0";
  $dir2 = opendir($imagedir);
  while( $file = readdir( $dir2 ) ) {
    if (eregi("runidl.sh", $file)) {
      $runidlexists = "1";
    }
  }
  if(! $runidlexists) {
    Exec("cd $imagedir;
          echo '#!/bin/sh' > runidl.sh;
          echo '' >> runidl.sh;
          echo 'IDL_PATH=../../../Idl:/usr/local/rsi/idl/lib:\${IDL_PATH}' >> runidl.sh;
          echo 'IDL_STARTUP=../../../Idl/idlrc_gui' >> runidl.sh;
          echo 'export IDL_PATH IDL_STARTUP' >> runidl.sh;
          echo '' >> runidl.sh;
          echo 'idl batch$macroextension' >> runidl.sh;
          chmod 755 runidl.sh");
  }
  $anyvectors = "";
  $uvectordir = "";
  if($uvector == "NO") {
    $uvector = "";
  } else {
    $anyvectors = "multiplot=1\n";
    if ($slice == "X") { $uvectordir = " uy;uz"; }
    if ($slice == "Y") { $uvectordir = " ux;uz"; }
    if ($slice == "Z") { $uvectordir = " ux;uy"; }
  }
  $bvectordir = "";
  if($bvector == "NO") {
    $bvector = "";
  } else {
    $anyvectors = "multiplot=1\n";
    if ($slice == "X") { $bvectordir = " by;bz"; }
    if ($slice == "Y") { $bvectordir = " bx;bz"; }
    if ($slice == "Z") { $bvectordir = " bx;by"; }
  }
  $jvectordir = "";
  if($jvector == "NO") {
    $jvector = "";
  } else {
    $anyvectors = "multiplot=1\n";
    if ($slice == "X") { $jvectordir = " jy;jz"; }
    if ($slice == "Y") { $jvectordir = " jx;jz"; }
    if ($slice == "Z") { $jvectordir = " jx;jy"; }
  }
  $funcvectors = "${uvectordir}${bvectordir}${jvectordir}";
  $modevectors = "${uvector}${bvector}${jvector}";
  if($view) {
    if(! $viewx1) { $viewx1 = "-1."; }
    if(! $viewx2) { $viewx2 = "1."; }
    if(! $viewy1) { $viewy1 = "-1."; }
    if(! $viewy2) { $viewy2 = "1."; }
	 $view = "!x.range = [$viewx1, $viewx2]\n!y.range = [$viewy1, $viewy2]\n";
  }
  $cwd = getcwd();
  $filename = "filename='file.out'\n";
  $func = "func='$variable1$funcvectors'\n";
  $plottitle = "plottitle='$variable1$funcvectors'\n";
  $plotmode = "plotmode='$conttype1$grid$modevectors'\n";
  $transform = "transform='n'\n";
  $range = "";
  if($contrange1 == "CUSTOM") {
    if(abs($contmin1) < 0.00001) { $contmin1 = "-0.00001"; }
    $range = "autorange='n'\nfmin=[$contmin1]\nfmax=[$contmax1]\n";
  }
  $saveps = "set_device,'image.ps',/port\n.r animate\nclose_device\n";
  Exec("mkdir $batchdir/$tmpdir;
        rsync -av $imagedir/runidl.sh $batchdir/$tmpdir/;
        cd $batchdir/$tmpdir;
        ln -s $cwd/$filedir/$plotfile file.out;
        echo \"$filename$view$func$plottitle$plotmode$anyvectors$transform$range$saveps\" > batch$macroextension");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot2() {   // Save defaults file
  global $imagedir, $number;

  $parameter = array('variable1', 'conttype1', 'contrange1', 'contmin1', 'contmax1', 'slice', 'grid', 'uvector', 'bvector', 'jvector', 'view', 'viewx1', 'viewx2', 'viewy1', 'viewy2');

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
  global $batchdir, $tmpdir, $file1;

  Exec("cd $batchdir/$tmpdir;
        echo '#!/bin/sh' > batchscript.sh;
        echo '' >> batchscript.sh;
        echo './runidl.sh >& batch.log' >> batchscript.sh;
        echo 'cp -f batch.log ../.' >> batchscript.sh;
        echo 'mv *ps $file1' >> batchscript.sh;
        chmod 755 batchscript.sh");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot4() {   // Download for self-plot instructions
  global $batchdir, $tmpdir, $filedir, $plotfile, $imagedir, $number, $macroextension;

  Exec("cd $batchdir/$tmpdir;
        echo '<h3>Download files to recreate this figure manually.</h3>' >> page.html;
        echo '<A HREF=\"${filedir}/${plotfile}\">${plotfile}</A><br>' >> page.html;
        echo '<A HREF=\"${imagedir}/batch-${number}${macroextension}\">batch-${number}${macroextension}</A><br>' >> page.html;
        echo '<A HREF=\"${imagedir}/runidl.sh\">runidl.sh</A><br>' >> page.html;
        echo '<br>' >> page.html;
        echo 'Then fix the path in runidl.sh and run<br>' >> page.html;
        echo '<CENTER>ln -s ${plotfile} file.out; ln -s batch-${number}${macroextension} batch${macroextension}; ./runidl.sh</CENTER>' >> page.html;
        echo '' >> page.html");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ?>
