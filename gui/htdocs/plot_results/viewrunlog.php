<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php include("../site_header_empty.php"); ?>

<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname', 'logfile');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html><head>
<title>SWMF GUI: Show runlog</title>
</head><body>
<div class="indent">
<pre>

<?php
  echo "$logfile for run: $runname\n\n\n";
  $return = "";
  Exec("cat ../plots/PLOT_$runname/$logfile", $return);
  foreach ($return as $tmp) {
    echo "$tmp\n";
  }
?>

</pre>
</div>
</body>
