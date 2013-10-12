<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php include("../site_header_empty.php"); ?>

<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname', 'refresh');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; }
?>

<?php
  echo "<html><head>";
  if($refresh) {
    echo "<META HTTP-EQUIV=\"refresh\" content=\"$refresh;URL=taillog.php?runname=$runname&refresh=$refresh\">";
  }
  $common = "taillog.php?runname=$runname&refresh=";
  echo "
</head><body>
&nbsp <b>View full <a href=\"viewlog.php?runname=$runname&logfile=runlog#end\" TARGET=\"_log\">
 runlog</a> file.</b><br>
&nbsp <b>Tail of runlog refresh rate:</b>
<select onchange=\"window.location='$common'+this.value\">
<option value=\"\">Never</option>
  ";
  $choices = array('2', '10', '60');
  foreach ($choices as $choice) {
    if($refresh==$choice) {
      echo "<option value=\"$choice\" SELECTED>$choice Seconds</option>";
    } else {
      echo "<option value=\"$choice\"         >$choice Seconds</option>";
    }
  }
  echo "
</select><br>
<div class=\"indent\">
<pre>\n";
  $return = "";
  Exec("tail -10 ../runs/RUN_$runname/runlog", $return);
  foreach ($return as $tmp) {
    echo "$tmp\n";
  }
?>
</pre></div></body>
