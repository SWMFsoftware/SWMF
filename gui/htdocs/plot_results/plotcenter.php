<?php

echo "<h1>SWMF GUI: Create Plots for \"$runname\"</h1><br>";

if (! $plotstyle) {  // Print initial info only
   PrintComponents();
   if ($cmp) {
     echo "<H2>Select component filetype from right column</H2>";
   } else {
     echo "<H2>Select component above</H2>";
     echo "<br><br><h3><a href=\"quicklook.php?runname=$runname\">Quick Look Plots</a></h3>";
     if( is_file("$runpath/runlog")) {
       echo "<br><br><h3><a href=\"viewrunlog.php?runname=$runname&logfile=runlog\"TARGET=\"_log\">View runlog</a></h3>";
     }
   }
   echo "<br>";
   echo '</BODY>';
   exit();
}


PrintComponents();

$cleanplotstyle = htmlentities(stripslashes($plotstyle), ENT_QUOTES);
$plotstyle = urldecode(stripslashes($plotstyle));

$plotfilelist = GetPlotList("$runpath/$cmp");
$countfiles = count($plotfilelist);
if($countfiles < 1) {
  echo "<H2>No plotfiles found!</H2><br>";
  echo '</BODY>';
  exit();
}

// If no plotfile is selected, choose the last available file.
if (! $plotfile) { SetLastPlotfile($plotfilelist, $plotfile); }

echo "<H2>Select desired plot options and 'Update Plot'</H2><br>";

// Load default and custom variables
include("$runpath/images/${cmp}_$plottype/defaultsBASE.php");

// Load plot specific form.
//
if($PlotApplication == "tecplot")    { include("plot_3Dplt.php"); }
if($PlotApplication == "tecplot1D")  { include("plot_1Dlog.php"); }
if($PlotApplication == "tecplot2D")  { include("plot_2Dplt.php"); }
if($PlotApplication == "tecplotLOS") { include("plot_LOS.php"); }
if($PlotApplication == "IEidl")      { include("plot_IEidl.php"); }
if($PlotApplication == "2Didl")      { include("plot_2Didl.php"); }
form1();

echo "</BODY>\n";

?>
