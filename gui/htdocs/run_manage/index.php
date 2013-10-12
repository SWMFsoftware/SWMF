<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Manage Runs</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Manage Runs</h1>
<BR CLEAR=ALL>

<?php
  if($runname) {
    echo "
<br>
<h2><a href=\"/run_manage/index.php\">Return</a> to manage runs without changes.</h2>
<br><br><br>
<h2>Run selected: $runname</h2>
<br>
<h2><a href=\"/run_manage/remove.php?runname=$runname\">REMOVE</a> (THERE IS NO UNDO ONCE CLICKED!)</h2>
<br><br><br>
    ";
    $return = "";
    $fname = "../runs/RUN_$runname/.codename";
    $fh = fopen($fname, "rt");
    $fcontent = fread($fh, filesize($fname));
    fclose($fh); 
    $codename = trim($fcontent);
    echo "
<h2>Run attached to code \"$codename\"</h2><br>
<h2>Run Directory Files</h2>
    ";
    $return = "";
    Exec("cd ../runs/RUN_$runname; ls -l > .files", $return);
    $fname = "../runs/RUN_$runname/.files";
    $fh = fopen($fname, "rt");
    $fcontent = fread($fh, filesize($fname));
    fclose($fh); 
    echo "
<pre>$fcontent</pre>
<BR CLEAR=ALL>
<hr>
    ";
  } else {
    echo "
<h2>View, Reuse, or Delete run from list:</h2><br>
<div class=\"indent\">
<table>
    ";
    $runs = "";
    $tmpdir = opendir( "$HOMEDIR/runs" );
    while( $run = readdir( $tmpdir ) ) {
      if (ereg("RUN_", $run)) { $runs[] = $run; }
    }
    if($runs){ sort($runs);
      foreach ($runs as $run) {
        $explode = explode("_", $run, 2); $name = trim($explode[1]);
        $inst = "Created"; $conf = ""; $comp = "";
        if (is_file("$HOMEDIR/runs/RUN2_${name}/core")) { $conf = "/ Setup & Run"; }
        if (is_file("$HOMEDIR/runs/RUN3_${name}/core")) { $conf = "/ Setup & Run"; $comp = "/ Completed"; }
        echo "
<tr>
<td><b>$name</b></td>
<td>&nbsp&nbsp</td>
<td>$inst</td><td>$conf</td><td>$comp</td>
<td>&nbsp&nbsp</td>
<td><a href=\"/run_manage/view.php?runname=$name\">View</a></td>
<td>&nbsp&nbsp</td>
<td><a href=\"/run_manage/reuse.php?runname=$name\">Reuse</a></td>
<td>&nbsp&nbsp</td>
<td><a href=\"/run_manage/index.php?runname=$name\">Delete</a></td>
</tr>
        ";
      }
    }
    echo "
</table>
</div>
    ";

  }
?>

<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>

</body>
</html>
