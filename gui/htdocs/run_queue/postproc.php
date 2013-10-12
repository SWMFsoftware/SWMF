<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php include('../config.php'); ?>

<?php    // Set posted value
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = $_POST[$name]; } ?>

<html>
<head>
<title>SWMF GUI: Postprocess Job</title>
<META HTTP-EQUIV=REFRESH CONTENT="0;URL=../plot_results/plot.php?runname=<?php echo $runname ?>">
</head>
<body>

<h3>This log window will clear shortly, please wait ...</h3>

<?php $pptime = `date +%Y%m%d%H%M%S`;
     $pptime = trim($pptime);
   Exec("cd ../runs;
         rm RUN2_$runname;
         ln -s RUN_$runname RUN3_$runname;
         cd RUN_$runname;
         PostProc.pl -o PLOT_${runname};
         rsync -a Idl PLOT_${runname}/;
         sleep 2");
   if(!  is_dir("$HOMEDIR/plots")) {
     Exec("mkdir $HOMEDIR/plots"); }
   Exec("cd $HOMEDIR/runs/RUN_$runname;
         mv PLOT_$runname ../../plots/.;
         ln -s ../../plots/PLOT_$runname plots_$pptime");
 ?>

</body>
