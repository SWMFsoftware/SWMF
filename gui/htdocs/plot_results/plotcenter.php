<?php

echo "<h1>SWMF GUI: Create Plots for \"$runname\"</h1><br>";

if (! $plotstyle) {  // Print initial info only
   PrintComponents();
   if ($cmp) {
     echo "<H2>Select component filetype from right column</H2>";
   } else {
     echo "<H2>Select component above</H2>";
     echo "<br><br><h3><a href=\"quicklook.php?runname=$runname&wait=1\"TARGET=\"_quick\">Quick Look Plots</a></h3>";
     if( is_file("$runpath/runlog")) {
       echo "<br><br><h3><a href=\"viewrunlog.php?runname=$runname&logfile=runlog\"TARGET=\"_log\">View runlog</a></h3>";
     }
   }
   echo "<br>";
   exit();
}


PrintComponents();

$cleanplotstyle = htmlentities(stripslashes($plotstyle), ENT_QUOTES);
$plotstyle = urldecode(stripslashes($plotstyle));

$plotfilelist = GetPlotList("$runpath/$cmp");
$countfiles = count($plotfilelist);
if($countfiles < 1) {
  echo "<H2>No plotfiles found!</H2><br>";
  exit();
}

// If no plotfile is selected, choose the last available file.
if (! $plotfile) { SetLastPlotfile($plotfilelist, $plotfile); }

echo "<H2>Select desired plot options and 'Update Plot'</H2><br>";

// Load default and custom variables
include("${runpath}/images/${cmp}_${plottype}/defaultsBASE.php");
include("plot_${loadfile}.php");
form1();

echo "<br><br>\n";

$plotslist = GetFileListByStyle("$runpath/images/${cmp}_${plottype}", "-${plotstyle}");
$countfiles = count($plotslist);
echo "<H2>Plots already created from selected stye: $countfiles</H2><br>";

$plotpath = "/plots/PLOT_$runname/images/${cmp}_${plottype}/";
if($countfiles > 0) {
  PrintPlotMenu($plotslist, "plotname", "$plotpath");
}

?>
