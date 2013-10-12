<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('cmp', 'runpath', 'plotfile');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Show pltview output</title>
</head>
<body>

<h3>pltview for currently selected file</h3>
<pre>
<?php
   $return = "";
   Exec("pltview $runpath/$cmp/$plotfile", $return);
      foreach ($return as $tmp) { echo "$tmp\n"; }
?>
</pre>

</body>
