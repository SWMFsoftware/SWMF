<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<META HTTP-EQUIV=REFRESH CONTENT="0;URL=../run_new/">
</head>
<body>

<?php
   set_time_limit(0);
   Exec("cd ../codes/CODE2_$codename;
         make rundir;
         mkdir runCLEAN;
         rsync -aL --exclude \"CVS\" run/ runCLEAN/;
         rsync -aL --exclude \"CVS\" GM/BATSRUS/Idl runCLEAN/;
         cd ..;
         rm CODE2_$codename;
         ln -s CODE_$codename CODE3_$codename");
   $deffile = "../codes/CODE_$codename/runCLEAN/job";
   if(!$fh = fopen($deffile, 'w')) {
     echo "<h1>The file $deffile cannot be opened</h1>";
   } else {
     fwrite($fh, "#!/bin/sh\n");
     fwrite($fh, "#PBS -S /bin/sh\n");
     fwrite($fh, "#PBS -l nodes=NN\n");
     fwrite($fh, "#PBS -q csem\n");
     fwrite($fh, "#PBS -N guijob\n");
     fwrite($fh, "\n");
     fwrite($fh, "# set environment variable\n");
     fwrite($fh, "export P4_SOCKBUFSIZE=200000\n");
     fwrite($fh, "\n");
     fwrite($fh, "# cd into the current directory, which should be the run directory\n");
     fwrite($fh, "cd \$PBS_O_WORKDIR\n");
     fwrite($fh, "\n");
     fwrite($fh, "# run code\n");
     fwrite($fh, "/usr/local/mpi/bin/mpirun -np NP -machinefile \$PBS_NODEFILE SWMF.exe >& runlog\n");
     fclose($fh);
   }
   $deffile = "../codes/CODE_$codename/runCLEAN/NN";
   if(!$fh = fopen($deffile, 'w')) {
     echo "<h1>The file $deffile cannot be opened</h1>";
   } else {
     fwrite($fh, "1");
     fclose($fh);
   }
   $deffile = "../codes/CODE_$codename/runCLEAN/NP";
   if(!$fh = fopen($deffile, 'w')) {
     echo "<h1>The file $deffile cannot be opened</h1>";
   } else {
     fwrite($fh, "2");
     fclose($fh);
   }
?>

</body>
