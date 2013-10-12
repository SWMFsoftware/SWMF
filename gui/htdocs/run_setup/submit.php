<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // Set posted value
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = $_POST[$name]; } ?>

<html>
<head>
<title>SWMF GUI: Submit Run Job</title>
<META HTTP-EQUIV=REFRESH CONTENT="1;URL=../run_queue/queue.php?runname=<?php echo $runname ?>">
</head>
<body>

<h3>This log window will clear shortly, please wait ...</h3>

<?php $subtime = `date +subtime%Y%m%d%H%M%S`;
   $fname = "../runs/RUN_$runname/NN";
   $fh = fopen($fname, "rt"); $NN = fread($fh, filesize($fname)); fclose($fh);
   $fname = "../runs/RUN_$runname/NP";
   $fh = fopen($fname, "rt"); $NP = fread($fh, filesize($fname)); fclose($fh);
   Exec("cd ../runs;
         rm RUN1_$runname;
         ln -s RUN_$runname RUN2_$runname;
         cd RUN_$runname;
         cat job | sed s/NN/'$NN'/ | sed s/NP/'$NP'/ > job2
         mv -f job2 job;
         qsub job >& $subtime");
 ?>

</body>
