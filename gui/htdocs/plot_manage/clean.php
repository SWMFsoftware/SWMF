<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('plotname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<META HTTP-EQUIV=REFRESH CONTENT="0;URL=index.php">
</head>
<body>

<?php
   $return = "";
   Exec("cd ../plots/PLOT_$plotname;
         rm -rf images", $return);
//   foreach ($return as $tmp) {
//      echo "$tmp\n";
//   }
?>

</body>
