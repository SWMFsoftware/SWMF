<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // Set posted value
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = $_POST[$name]; } ?>

<html>
<head>
<title>SWMF GUI: Run Job</title>
<META HTTP-EQUIV=REFRESH CONTENT="3;URL=../run_queue/queue.php?runname=<?php echo $runname ?>">
</head>
<body>

<h3>The run has started, opening run monitor page ...</h3>

<?php
   set_time_limit(0);
   $fname = "../runs/RUN_$runname/NP";
   $fh = fopen($fname, "rt");
   $NP = fread($fh, filesize($fname));
   fclose($fh);
   $subtime = `date +runtime%Y%m%d%H%M%S`;
   $command = "mpirun -np $NP SWMF.exe";
   Exec("cd ../runs;
	 rm RUN1_$runname;
	 ln -s RUN_$runname RUN2_$runname;
	 cd RUN_$runname;
	 echo $command >> $subtime");
   Exec("cd ../runs/RUN_$runname;
	 $command >> runlog &");
 ?>

</body>
