<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Setup Run</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Setup Run "<?php echo "$runname" ?>"</h1>
<BR CLEAR=ALL>

<?php
   $fname = "../runs/RUN_$runname/.codename";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh);
   $codename = trim($fcontent);
?>
<h2>Run attached to code "<?php echo $codename ?>"</h2><br>
<?php
   if(! is_dir("../codes/CODE_$codename")) {
     echo "
<b><font color=red>
ERROR! Code attatched to this run cannot be found.
Because of this, we cannot properly check the validity of the input files
and cannot allow this run to proceed.
</font></b>
     ";
	  exit;
   }
?>
<h2>Run Directory Files</h2>
<div class="indent">
<?php $return = "";
   Exec("cd ../runs/RUN_$runname; ls -l > .files", $return);
   $fname = "../runs/RUN_$runname/.files";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh); 
?>
<pre><?php echo $fcontent ?></pre>
<br>
<FORM TARGET="" METHOD="post" ACTION="upload.php" ENCTYPE="multipart/form-data">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname ?>">
  <INPUT TYPE=hidden name=MAX_FILE_SIZE VALUE=100000>
  <INPUT TYPE=submit VALUE="Upload"> a new file to run directory (100kb limit) &nbsp&nbsp&nbsp
  <INPUT TYPE=file NAME=uploadedfile><br>
</FORM>
</div>
<BR CLEAR=ALL>
<hr>

<h2>Run Job</h2>
<br>
<div class="indent">
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Verify number of nodes/processors to use.</h3>
<?php
   $fname = "../runs/RUN_$runname/NN";
   $fh = fopen($fname, "rt"); $NN = fread($fh, filesize($fname)); fclose($fh);
   $fname = "../runs/RUN_$runname/NP";
   $fh = fopen($fname, "rt"); $NP = fread($fh, filesize($fname)); fclose($fh);
?>
<FORM TARGET="" METHOD="post" ACTION="savefile.php">
  <INPUT TYPE=hidden name=filename value="NN">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname ?>">
  Using <b><?php echo $NN ?></b> node(s).  Change value to <INPUT type=text size=4 name=data value="<?php echo $NN ?>">
  <INPUT TYPE=submit VALUE="save NN"></FORM>
<FORM TARGET="" METHOD="post" ACTION="savefile.php">
  <INPUT TYPE=hidden name=filename value="NP">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname ?>">
  Using <b><?php echo $NP ?></b> processor(s).  Change value to <INPUT type=text size=4 name=data value="<?php echo $NP ?>">
  <INPUT TYPE=submit VALUE="save NP"></FORM>
<?php
  if($NP < 1) { 
    echo "<pre><b><font color=red> ERROR, you need to use a positive integer number of processors!</font></b></pre>";
  }
?>

<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Input files need to pass parameter check before running.</h3>
<?php 
   Exec("cd ../codes/CODE_$codename;
         Scripts/TestParam.pl -n=$NP ../../runs/RUN_$runname/PARAM.in >& ../../runs/RUN_$runname/.testparam");
   $fname = "../runs/RUN_$runname/.testparam";
   $fcontent = "";
   if ( filesize($fname) > 0 ) {
     $fh = fopen($fname, "rt");
     $fcontent = fread($fh, filesize($fname));
     fclose($fh);
   }
   if (! $fcontent == "") {
	  echo "
<pre>
TestParam.pl output:
<b><font color=red> $fcontent </font></b>
</pre>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp You need to fix TestParam.pl errors before you can run.</h3>
     ";
   } else {
     echo "
<pre>
TestParam.pl output:
<b><font color=green> No errors were found. </font></b>
</pre>
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Make sure that the files below are correct before starting job!</h3>
<br>
<FORM TARGET=\"\" METHOD=\"post\" ACTION=\"submit.php\">
  <INPUT TYPE=hidden name=runname value=\"$runname\">
  <INPUT TYPE=submit VALUE=\"QSUB\"><b>&nbsp&nbsp Selecting QSUB will execute 'qsub job'</b></FORM>
<FORM TARGET=\"\" METHOD=\"post\" ACTION=\"run.php\">
  <INPUT TYPE=hidden name=runname value=\"$runname\">
  <INPUT TYPE=submit VALUE=\"EXEC\"><b>&nbsp&nbsp&nbsp Selecting EXEC will execute 'mpirun -np $NP SWMF.exe >& runlog'</b></FORM>
     ";
   }
?>
</div>
<BR CLEAR=ALL>

<!--
<hr>
<h2>Edit LAYOUT.in</h2>
<div class="indent">
<FORM TARGET="" METHOD="post" ACTION="saveLAYOUT.php">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname; ?>">
  <TABLE ALIGN="LEFT" BORDER=0 CELLPADDING=2 CELLSPACING=0>
  <tr><td>Name</td><td>First</td><td>Last</td><td>Stride</td></tr>
<?php
  $return = "";
  Exec("cat ../runs/RUN_$runname/LAYOUT.in", $return);
  $i = "0";
  $readit = "";
  foreach ($return as $tmp) {
    if(eregi("COMPONENTMAP", $tmp)) { $readit = "1";
	 } elseif(eregi("END", $tmp)) { $readit = "";
	 } elseif($readit) {
      $i += 1;
      $chunks = split(" ", $tmp, 2); $name   = trim($chunks[0]); $tmp = trim($chunks[1]);
      $chunks = split(" ", $tmp, 2); $first  = trim($chunks[0]); $tmp = trim($chunks[1]);
      $chunks = split(" ", $tmp, 2); $last   = trim($chunks[0]); $tmp = trim($chunks[1]);
      $chunks = split(" ", $tmp, 2); $stride = trim($chunks[0]);
      if($i > 1) { echo "</tr>"; }
      echo " 
  <INPUT TYPE=hidden name=name$i value=\"$name\"><TR>
  <td>$name</td>
  <td><INPUT type=text size=5 name=first$i value=\"$first\"></td>
  <td><INPUT type=text size=5 name=last$i value=\"$last\"></td>
  <td><INPUT type=text size=5 name=stride$i value=\"$stride\"></td>
      ";
	 }
  }
  echo "  <td><INPUT TYPE=submit VALUE=\"save LAYOUT.in\"></td></tr></TABLE></FORM>\n";
?>
</div>
<BR CLEAR=ALL><br>
--!>

<hr>
<h2>Edit LAYOUT.in</h2>
<br>
<div class="indent">
<?php $fname = "../runs/RUN_$runname/LAYOUT.in";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh) ?>
<FORM TARGET="" METHOD="post" ACTION="savefile.php">
  <INPUT TYPE=hidden name=filename value="LAYOUT.in">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname ?>">
  <textarea name=data rows=8 cols=82><?php echo $fcontent ?></textarea><br>
  <INPUT TYPE=submit VALUE="save LAYOUT.in"></FORM>
</div>
<BR CLEAR=ALL>

<hr>
<h2>Edit PARAM.in</h2>
<br>
<div class="indent">
<?php $fname = "../runs/RUN_$runname/PARAM.in";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh) ?>
<FORM TARGET="" METHOD="post" ACTION="savefile.php">
  <INPUT TYPE=hidden name=filename value="PARAM.in">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname ?>">
  <textarea name=data rows=30 cols=82><?php echo $fcontent ?></textarea><br>
  <INPUT TYPE=submit VALUE="save PARAM.in"></FORM>
</div>
<BR CLEAR=ALL>

<hr>
<h2>Edit job submission file</h2>
<br>
<div class="indent">
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp NOTE: NN and NP will be replaced with <?php echo $NN ?> and <?php echo $NP ?> when the job is submitted.</h3>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Be sure to pipe run output into the file "runlog"</h3>
<?php $fname = "../runs/RUN_$runname/job";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh) ?>
<FORM TARGET="" METHOD="post" ACTION="savefile.php">
  <INPUT TYPE=hidden name=filename value="job">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname ?>">
  <textarea name=data rows=15 cols=82><?php echo $fcontent ?></textarea><br>
  <INPUT TYPE=submit VALUE="save job"></FORM>
</div>
<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>

</body>
</html>
