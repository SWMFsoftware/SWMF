<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Config.pl help</title>
</head>
<body>

<pre>
<?php
   echo "Config.pl -h\n\n";
   $return = "";
   Exec("cd ../codes/CODE1_$codename; ./Config.pl -h", $return);
   foreach ($return as $tmp) {
      echo "$tmp\n";
   }
?>
</pre>

</body>
