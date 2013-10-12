<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php include("../site_header_empty.php"); ?>

<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename', 'refresh');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; }
?>

<?php
  echo "<html><head>";
  if($refresh) {
    echo "<META HTTP-EQUIV=\"refresh\" content=\"$refresh;URL=taillog.php?codename=$codename&refresh=$refresh\">";
  }
  $common = "taillog.php?codename=$codename&refresh=";
  echo "
</head><body>
&nbsp <b>View full <a href=\"viewlog.php?codename=$codename&logfile=compile.log#end\" TARGET=\"_log\">
 compile log</a> file.</b><br>
&nbsp <b>Tail of compile log refresh rate:</b>
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
  if(is_file("../codes/CODE_$codename/compile.log")) {
    $return = "";
    Exec("tail -10 ../codes/CODE_$codename/compile.log", $return);
    foreach ($return as $tmp) {
      echo "$tmp\n";
    }
  }
?>
</pre></div></body>
