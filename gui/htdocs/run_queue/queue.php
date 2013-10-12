<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Queue Run</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Queue Run "<?php echo "$runname" ?>"</h1>
<BR CLEAR=ALL>

<?php $return = "";
   $fname = "../runs/RUN_$runname/.codename";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh); 
?>
<h2>Run attached to code "<?php echo $fcontent ?>"</h2><br>

<h2>Run Directory Files</h2>
<?php $return = "";
   Exec("cd ../runs/RUN_$runname; ls -l > .files", $return);
   $fname = "../runs/RUN_$runname/.files";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh); 
?>
<pre><?php echo $fcontent ?></pre>
<BR CLEAR=ALL>
<hr>

<h2>Monitor Jobs</h2>
<?php
  $submits = array();
  $tmpdir = opendir( "../runs/RUN_$runname" );
  while( $submit = readdir( $tmpdir ) ) {
     if (ereg("subtime", $submit)) { $submits[] = $submit; }
  }
  if($submits){ sort($submits);
     foreach ($submits as $submit) {
        echo "
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Contents of $submit:</h3>
	";
	$fname = "../runs/RUN_$runname/$submit";
   	$fh = fopen($fname, "rt");
   	$fcontent = fread($fh, filesize($fname));
   	fclose($fh);
	echo "<pre>";
	echo $fcontent;
	echo "</pre>";
     }
     if (! is_file("../runs/RUN_$runname/SWMF.SUCCESS")) {
       echo "
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Output from 'qstat -a'</h3>
<pre>
       ";
       $return = "";
       Exec("qstat -a", $return);
       foreach ($return as $tmp) {
         echo "$tmp<br>";
       }
       echo "
<br>
</pre>
       ";
     }
  }
  $runs = array();
  $tmpdir = opendir( "../runs/RUN_$runname" );
  while( $run = readdir( $tmpdir ) ) {
     if (ereg("runtime", $run)) { $runs[] = $run; }
  }
  if($runs){ sort($runs);
     foreach ($runs as $run) {
        echo "
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Contents of $run:</h3>
	";
	$fname = "../runs/RUN_$runname/$run";
   	$fh = fopen($fname, "rt");
   	$fcontent = fread($fh, filesize($fname));
   	fclose($fh);
	echo "<pre>";
	echo $fcontent;
	echo "</pre>";
     }
     if (! is_file("../runs/RUN_$runname/SWMF.SUCCESS")) {
       echo "
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Output from 'ps -uxc | grep SWMF.exe | grep -v grep'</h3>
<pre>";
       $return = "";
       Exec("ps -uxc | grep SWMF.exe | grep -v grep", $return);
       echo "<br>";
       foreach ($return as $tmp) {
         echo "$tmp<br>";
       }
       echo "
<br>
</pre>
       ";
     }
  }
  if (is_file("../runs/RUN_$runname/SWMF.SUCCESS")) {
    echo "
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Job has finished.</h3>
<br>
    ";
  }
?>
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Code output runlog.</h3>
<br>
<?php
  if (! is_file("../runs/RUN_$runname/runlog")) {
    echo "
Job has not started.  No runlog exists.
    ";
  } else {
    echo "
<iframe width=105% height=210 frameborder=1 src=\"taillog.php?runname=$runname\"></iframe>
<br><br>
    ";
  }
?>
<br><br>
<hr>

<h2>Postprocess Files and Prepare To Visualize</h2>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Files moved to plots_[timestamp]</h3>
<br>
<?php
  if (! is_file("../runs/RUN_$runname/SWMF.SUCCESS")) {
    echo "
Job has not completed or has died.  Check runlog.
<br><br>
    ";
  } else {
    echo "
Job is done.  You may now PostProcess output.
<br><br>
<FORM TARGET=\"\" METHOD=\"post\" ACTION=\"postproc.php\">
  <INPUT TYPE=hidden name=runname value=\"$runname\">
  <INPUT TYPE=submit VALUE=\"PostProcesss\"></FORM>
    ";
  }
?>
<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>

</body>
</html>
