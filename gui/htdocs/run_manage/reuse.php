<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Reuse Run</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Reuse Run "<?php echo "$runname" ?>"</h1>
<BR CLEAR=ALL>

<br>
<h2><a href="/run_manage/index.php">Return</a> to manage runs without changes.</h2>
<br><br><br>

<?php $return = "";
   $fname = "../runs/RUN_$runname/.codename";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh); 
   $codename = trim($fcontent);
?>
<h2>Run attached to code "<?php echo $codename ?>"</h2><br>

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

<h2>Select This Job To Reuse</h2>
<br>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Create a unique run descriptor (no spaces) to describe your configuration, like Bz_flip.</h3>
<br>
<FORM TARGET="" ACTION="reuserundescriptor.php">
  <INPUT TYPE=hidden name=codename value="<?php echo $codename; ?>">
  <INPUT TYPE=hidden name=runname value="<?php echo $runname; ?>">
  <LABEL>Run Descriptor: </LABEL><INPUT type=text name="descriptor">
  <INPUT TYPE=submit VALUE="Create New Run"> attached to code "<?php echo $codename; ?>"</FORM>
<BR CLEAR=ALL>
<hr>

<h2>Nodes & Processors</h2>
<br>
<?php $fname = "../runs/RUN_$runname/NN";
   $fh = fopen($fname, "rt");
   $NN = fread($fh, filesize($fname));
   fclose($fh) ?>
<?php $fname = "../runs/RUN_$runname/NP";
   $fh = fopen($fname, "rt");
   $NP = fread($fh, filesize($fname));
   fclose($fh) ?>
  Nodes (NN) = <?php echo $NN ?>, Processors (NP) = <?php echo $NP ?><br>
<BR CLEAR=ALL>

<h2>Job submission file</h2>
<br>
<?php $fname = "../runs/RUN_$runname/job";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh) ?>
  <textarea name=data rows=15 cols=82><?php echo $fcontent ?></textarea><br>
<BR CLEAR=ALL>

<h2>LAYOUT.in</h2>
<br>
<?php $fname = "../runs/RUN_$runname/LAYOUT.in";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh); ?>
  <textarea name=data rows=8 cols=82><?php echo $fcontent ?></textarea><br>
<BR CLEAR=ALL>

<h2>PARAM.in</h2>
<br>
<?php $fname = "../runs/RUN_$runname/PARAM.in";
   $fh = fopen($fname, "rt");
   $fcontent = fread($fh, filesize($fname));
   fclose($fh) ?>
  <textarea name=data rows=30 cols=82><?php echo $fcontent ?></textarea><br>
<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>

</body>
</html>
