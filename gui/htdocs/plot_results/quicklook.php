<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php
  $parameter = array('runname', 'wait');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; }

  if (($wait)) {
    echo "
<html>
<head>
<title>GUI Quicklook: ... rendering ...</title>
<META HTTP-EQUIV=\"refresh\" content=\"0;URL=quicklook.php?runname=${runname}\">
</head>
<body>
<br><center><img src=\"../images/pleasewait.gif\"></center>
</body>
</html>
    ";
    exit();
  }

  echo "
<html>
<head>
<title>GUI Quicklook</title>
</head>
<body>
<h2>Quick Look plots for run $runname</h2>
<table cellpadding=\"5\">
  ";

  $column = "0";
  include("plotfunctions.php");

  $runpath = "../plots/PLOT_$runname";

  $dir = opendir( "$runpath" );
  while( $file = readdir( $dir ) ) {
    if ( ereg("[A-Z][A-Z]", $file) && is_dir("$runpath/$file") ) {
      $cmp = $file;
      $filedir = "$runpath/$cmp";

      if(! is_dir("$runpath/images")) {
        Exec("cd $runpath; mkdir images");
      }
      Exec("rsync -a --exclude \"CVS\" BASE/$cmp* $runpath/images/");
      $dir2 = opendir( "$runpath/images" );
      while( $file = readdir( $dir2 ) ) {
        if (ereg("$cmp", $file)) {
          $pieces = explode("_", $file);
          $plottype = $pieces[1];

          include("$runpath/images/${cmp}_$plottype/defaultsBASE.php");

          $number = "001";

          $plotfile = "";
          $plotfilelist = GetPlotList("$runpath/$cmp");
          $countfiles = count($plotfilelist);
          if($countfiles > 0) {
            if(! $plotfile) { SetLastPlotfile($plotfilelist, $plotfile); }
            include("makeplot.php");
            $column++;
            if("$column" == "3") { $column = "1"; }
            if("$column" == "1") { echo "<tr>"; }
            echo "
<td>
<center>
<b>${cmp} &nbsp&nbsp ${plottype}: ${plotfileclip}</b><br><br>
<IMG SRC=\"$imagedir/$file2\" width=95% BORDER=0>
</center>
</td>
            ";
            if("$column" == "2") { echo "</tr>"; }
          }

        }
      }

    }
  }
  if("$column" == "1") {
    echo "</tr>";
  }
?>
</table>
</body>
</html>
