<?php

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function form1() {   //
  global $runpath, $runname, $cmp, $plottype, $plotstyle, $plotfile, $plotfilelist;
  global $countfiles;

  include("$runpath/images/${cmp}_$plottype/defaults-$plotstyle.php");

  echo "
<b>Plotfile:</b> T=Year:Month:Day - Hour:Min:Sec ($countfiles files found)<br>
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

  // CONTOUR
  echo "
<b>Contour</b><br>
<div class=\"indent\">
<table>
  ";
  $numvars = count($variables);
  if ($numvars > 6) {
    echo "<b>NOTE: Not expecting more than 6 variables in file, fix form</b>";
//    exit();
  }
  $i = 0;
  foreach ($variables as $value) {
    $i++;
    if ($i <= 2) {
      echo "
<INPUT TYPE=hidden name=contplot${i} value=\"N\">
<INPUT TYPE=hidden name=contrange${i} value=\"MM\">
<INPUT TYPE=hidden name=contrangevalue${i} value=\"\">
      ";
      //    } else {
    } else if ($i <= 6) {
      // Set CHECKED for defcont
      $defcontY = ""; $defcontN = " CHECKED";
      $tmpname = "defcontplot$i";
      if($$tmpname == "Y") {$defcontY = " CHECKED"; $defcontN = "";}
      // Set CHECKED for defcontrange
      $defcontrangeMM = " CHECKED"; $defcontrangeCUSTOM = "";
      $tmpname = "defcontrange$i";
      if($$tmpname == "CUSTOM") {$defcontrangeMM = ""; $defcontrangeCUSTOM = " CHECKED";}
      // Set value for defcontrangevalue
      $tmpname = "defcontrangevalue$i";
      $tmpvalue = $$tmpname;
      echo "
<tr>
<td>$value</td>
<td>&nbsp&nbsp</td>
<td><INPUT TYPE=radio name=contplot${i} value=\"Y\"$defcontY>Yes</td>
<td><INPUT TYPE=radio name=contplot${i} value=\"N\"$defcontN>No</td>
<td>&nbsp&nbsp</td>
<td>Range:</td>
<td><INPUT TYPE=radio name=contrange${i} value=\"MM\"$defcontrangeMM>Max</td>
<td><INPUT TYPE=radio name=contrange${i} value=\"CUSTOM\"$defcontrangeCUSTOM>Custom</td>
<td><INPUT type=text size=5 name=contrangevalue${i} value=\"$tmpvalue\"></td>
</tr>
      ";
    }
  }
  echo "
</table>
</div>
 <br>
<INPUT TYPE=submit VALUE=\"Update Plot\"></FORM>
  ";
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot0() {   // Check for errors
  $parameter = array('contplot1', 'contplot2', 'contplot3', 'contplot4', 'contplot5', 'contplot6');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  foreach ($parameter as $name) {
    if ($$name == "Y") { return 0; }
  }
  return 1;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot1() {   // Create scripts
  global $imagedir, $filedir, $plotfile, $batchdir, $tmpdir, $macroextension;

  $parameter = array('contplot1', 'contplot2', 'contplot3', 'contplot4', 'contplot5', 'contplot6', 'contrange1', 'contrange2', 'contrange3', 'contrange4', 'contrange5', 'contrange6', 'contrangevalue1', 'contrangevalue2', 'contrangevalue3', 'contrangevalue4', 'contrangevalue5', 'contrangevalue6');
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
          echo 'idl < batch$macroextension' >> runidl.sh;
          chmod 755 runidl.sh");
  }
  $i = 0;
  $stringvar = "";
  $stringrange = "";
  foreach ($parameter as $name) {
    if(eregi("contplot", $name)) {
      $j = $i;
      $i++;
      $tmpname = "contrange$i";
      $tmpname2 = "contrangevalue$i";
      if($$tmpname == "MM") {$$tmpname2 = "";}
      if ($$name == "Y") {
         $stringvar = "$stringvar$j\n";
	 $range = $$tmpname2;
         $stringrange = "$stringrange$range\n";
      }
    }
  }
  $cwd = getcwd();
  Exec("mkdir $batchdir/$tmpdir;
        rsync -av $imagedir/runidl.sh $batchdir/$tmpdir/;
        cd $batchdir/$tmpdir;
        ln -s $cwd/$filedir/$plotfile file.idl;
        echo \".r plot_iono_gui\nfile.idl\n\n\n\n$stringvar\n$stringrange\" > batch$macroextension");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function plot2() {   // Save defaults file
  global $imagedir, $number;

  $parameter = array('contplot1', 'contplot2', 'contplot3', 'contplot4', 'contplot5', 'contplot6', 'contrange1', 'contrange2', 'contrange3', 'contrange4', 'contrange5', 'contrange6', 'contrangevalue1', 'contrangevalue2', 'contrangevalue3', 'contrangevalue4', 'contrangevalue5', 'contrangevalue6');

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
