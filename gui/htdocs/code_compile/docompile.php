<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename', 'options');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<META HTTP-EQUIV=REFRESH CONTENT="5;URL=compile.php?codename=<?php echo $codename ?>">
</head>
<body>

<h2>The code is being compiled in the background.</h2>
<h2>This page will refresh to the compile options page in 5 seconds</h2>
<h2>Follow 'compile log' link to check on compilation progress.</h2>

<?php
   set_time_limit(0);
   $command = "make $options";
   $return = "";
   Exec("cd ../codes/CODE2_$codename;
         echo $command >> command.log;
         echo \" \" >> compile.log;
         echo \"===NEXT COMMAND===\" >> compile.log;
         echo $command >> compile.log;
         echo \" \" >> compile.log;
         $command >> compile.log &", $return);
//   foreach ($return as $tmp) {
//      echo "$tmp\n";
//   }
?>

</body>
