<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php include("../site_header_empty.php"); ?>

<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname', 'refresh', 'logfile');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<?php
  echo "<html><head>";
  echo "<title>SWMF GUI: Show logs</title>";
  if($refresh) {
    echo "<META HTTP-EQUIV=\"refresh\" content=\"$refresh;URL=viewlog.php?runname=$runname&refresh=$refresh&logfile=$logfile#end\">";
  }
  echo "</head><body>";
  echo "<div class=\"indent\">";
  echo "<pre>";
  echo "cat $logfile\n\n";
  $return = "";
  Exec("cat ../runs/RUN_$runname/$logfile", $return);
  foreach ($return as $tmp) {
    echo "$tmp\n";
  }
?>
</pre>
</div>
<div id="end"></div><br>
<?php
  $common = "viewlog.php?runname=$runname&logfile=$logfile&refresh=";
  echo "
&nbsp <b>Log refresh rate:</b>
<select onchange=\"window.location='$common'+this.value\">
<option value=\"#end\">Never</option>
  ";
  $choices = array('2', '10', '60');
  foreach ($choices as $choice) {
    if($refresh==$choice) {
      echo "<option value=\"$choice#end\" SELECTED>$choice Seconds</option>";
    } else {
      echo "<option value=\"$choice#end\"         >$choice Seconds</option>";
    }
  }
  echo "</select><br>";
?>
</body>
