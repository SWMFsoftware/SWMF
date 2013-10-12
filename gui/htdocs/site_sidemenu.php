<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php echo "<!-- BEGIN " . basename(__FILE__) . " -->\n"; ?>

<div id="leftnavbar">
<div class="navlist">

<?php 
   $DOcode = "1";
   $DOrun = "1";
   if(is_file("$HOMEDIR/../../Config.pl")) { $IScode = "1"; } else { $IScode = "0"; }
   if(!$IScode) { $DOcode = "0"; $DOrun = "0"; }
   $dir = getcwd();
   echo "
<br>
<ul>
  <h3>SWMF GUI</h3>
</ul>
<ul>
  <li><a href=\"/\"title=\"home\">Home</a></li>
   ";
   if ($IScode) {
     if(is_dir("$HOMEDIR/Manuals")) {
       echo "
  <li><a href=\"/Manuals/SWMF/\"title=\"online SWMF manual\">SWMF Manual</a></li>
  <li><a href=\"/Manuals/REFERENCE/\"title=\"online REFERENCE manual\">REFERENCE Manual</a></li>
       ";
     } else {
       echo "
  No Manuals Created
       ";
     }
   }
   if ($DOcode) {
     echo "
  <hr>
  <li><a href=\"/code_install/\"title=\"create new code directory\">Create Code Directory</a></li>
  <li><a href=\"/code_config/\"title=\"configure created code\">Configure Code</a></li>
     ";
     if (eregi("code_config", $dir)) {
		  $list = FindList("$HOMEDIR/codes", "CODE1");
		  echo "<ul>\n";
		  foreach ($list as $item) { echo "<li><a href=\"/code_config/configure.php?codename=$item\">$item</a></li>\n"; }
		  echo "     </ul>\n";
     }
     echo "
  <li><a href=\"/code_compile/\"title=\"compile configured code\">Compile Code</a></li>
     ";
     if (eregi("code_compile", $dir)) {
		  $list = FindList("$HOMEDIR/codes", "CODE2");
		  echo "<ul>\n";
		  foreach ($list as $item) { echo "<li><a href=\"/code_compile/compile.php?codename=$item\">$item</a></li>\n"; }
		  echo "     </ul>\n";
     }
     echo "
  <li><a href=\"/code_manage/\"title=\"manage installed codes: delete, etc.\">Manage Codes</a></li>
     ";
   }
   if ($DOrun) {
     echo "
  <hr>
  <li><a href=\"/run_new/\"title=\"create new run directory and assign to code\">Create Run Directory</a></li>
  <li><a href=\"/run_setup/\"title=\"setup input files that define run and execute\">Setup & Execute Run</a></li>
     ";
     if (eregi("run_setup", $dir)) {
		  $list = FindList("$HOMEDIR/runs", "RUN1");
		  echo "<ul>\n";
		  foreach ($list as $item) { echo "<li><a href=\"/run_setup/setup.php?runname=$item\">$item</a></li>\n"; }
		  echo "     </ul>\n";
     }
     echo "
  <li><a href=\"/run_queue/\"title=\"monitor running jobs and postprocess when finished\">Monitor & Process Run</a></li>
     ";
     if (eregi("run_queue", $dir)) {
		  $list = FindList("$HOMEDIR/runs", "RUN2");
		  echo "<ul>\n";
		  foreach ($list as $item) { echo "<li><a href=\"/run_queue/queue.php?runname=$item\">$item</a></li>\n"; }
		  echo "     </ul>\n";
     }
     echo "
  <li><a href=\"/run_manage/\"title=\"manage runs: delete, reuse, etc.\">Manage Runs</a></li>
     ";
  }
  echo "
  <hr>
  <li><a href=\"/plot_results/\"title=\"create plots from runs and data\">Create Plots</a></li>
  ";
  if (eregi("plot_results", $dir)) {
	  $list = FindList("$HOMEDIR/plots", "PLOT");
	  echo "<ul>\n";
	  foreach ($list as $item) { echo "<li><a href=\"/plot_results/plot.php?runname=$item\">$item</a></li>\n"; }
	  echo "     </ul>\n";

	  $list = FindList("$HOMEDIR/plots", "DATA");
	  echo "<ul>\n";
	  foreach ($list as $item) { echo "<li><a href=\"/plot_results/dataplot.php?runname=$item\">$item</a></li>\n"; }
	  echo "     </ul>\n";
  }
  echo "
  <li><a href=\"/plot_manage/\"title=\"manage installed plots: delete, clean, etc.\">Manage Plots</a></li>
  ";
  echo "
  <hr>
<!--   <li><a href=\"/Perry/\"title=\"Slice-N-Dice\">Perry's Slices</a></li> -->
<!--   <hr> -->
  &nbsp
</ul>
<ul>
  <h3>LINKS</h3>
</ul>
<ul>
  <li><a href=\"http://csem.engin.umich.edu\" title=\"Comprehensive Space Environment Modeling\">CSEM</a></li>
</ul>
<center>
<br><br>
<img src=\"/images/SWMF.png\" width=\"160\">
<br><br>
<img src=\"/images/CSEM.png\" width=\"100\">
<br><br>
<img src=\"/images/UM.gif\" width=\"150\">
</center>
</div>
</div>
  ";
 ?>

<?php echo "<!-- END " . basename(__FILE__) . " -->\n"; ?>

<?php
function FindList($where='.', $what) {
  $list = array();
  $tmpdir = opendir( "$where" );
  while( $item = readdir( $tmpdir ) ) {
    if (ereg("${what}_", $item)) { $explode = explode("_", $item, 2); $list[] = $explode[1]; }
  }
  sort($list);
  return $list;
}
?>
