<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename', 'install', 'options', 'GM', 'IE', 'IH', 'IM', 'PS', 'PW', 'RB', 'SC', 'SP', 'UA', 'compiler', 'mpiversion', 'grid', 'grid1', 'grid2', 'grid3', 'grid4', 'grid5', 'grid6', 'grid7', 'grid8', 'grid9');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<?php
  $vopt = "";
  if($GM) { if($vopt) { $vopt .= ",$GM"; } else { $vopt = "-v=$GM"; } }
  if($IE) { if($vopt) { $vopt .= ",$IE"; } else { $vopt = "-v=$IE"; } }
  if($IH) { if($vopt) { $vopt .= ",$IH"; } else { $vopt = "-v=$IH"; } }
  if($IM) { if($vopt) { $vopt .= ",$IM"; } else { $vopt = "-v=$IM"; } }
  if($PS) { if($vopt) { $vopt .= ",$PS"; } else { $vopt = "-v=$PS"; } }
  if($PW) { if($vopt) { $vopt .= ",$PW"; } else { $vopt = "-v=$PW"; } }
  if($RB) { if($vopt) { $vopt .= ",$RB"; } else { $vopt = "-v=$RB"; } }
  if($SC) { if($vopt) { $vopt .= ",$SC"; } else { $vopt = "-v=$SC"; } }
  if($SP) { if($vopt) { $vopt .= ",$SP"; } else { $vopt = "-v=$SP"; } }
  if($UA) { if($vopt) { $vopt .= ",$UA"; } else { $vopt = "-v=$UA"; } }
  $copt = ""; if($compiler) { $copt = "-compiler=$compiler"; }
  $mopt = ""; if($mpiversion) { $mopt = "-mpi=$mpiversion"; }
  $gopt = ""; if($grid) {
    $gopt = "-g={$grid}:$grid1";
    if($grid2) { $gopt .= ",$grid2"; }
    if($grid3) { $gopt .= ",$grid3"; }
    if($grid4) { $gopt .= ",$grid4"; }
    if($grid5) { $gopt .= ",$grid5"; }
    if($grid6) { $gopt .= ",$grid6"; }
    if($grid7) { $gopt .= ",$grid7"; }
    if($grid8) { $gopt .= ",$grid8"; }
    if($grid9) { $gopt .= ",$grid9"; }
  }
?>

<html>
<head>
<META HTTP-EQUIV=REFRESH CONTENT="0;URL=configure.php?codename=<?php echo $codename ?>">
</head>
<body>

<?php
   set_time_limit(0);
   $command = "./Config.pl {$install} {$copt} {$mopt} {$vopt} {$gopt} {$options}";
   $return = "";
   Exec("cd ../codes/CODE1_$codename;
         echo $command >> command.log;
         $command >& .tmp", $return);
   //print_r($return);
?>

</body>
