<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename', 'descriptor', 'runname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Reuse Run Directory</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Reuse Run Directory</h1>
<BR CLEAR=ALL>

<?php if($descriptor) {
     if($codename) {
        // Check if run already exists
        if (is_dir("../runs/RUN_$descriptor")) {
	  echo "<h2>Descriptor (run) already exists, try again.</h2><BR CLEAR=ALL>";
          include("../site_footer.php");
          echo "</body></html>";
          exit();
	}

        // Check if plot already exists
        if (is_dir("../plots/PLOT_$descriptor")) {
	  echo "<h2>Descriptor (plot) already exists, try again.</h2><BR CLEAR=ALL>";
          include("../site_footer.php");
          echo "</body></html>";
          exit();
	}

        $return = "";
        Exec("cd ../runs;
              mkdir RUN_$descriptor;
              rsync -av --exclude \"plots_*\" --exclude \"SWMF.SUCCESS\" --exclude \"subtime*\" --exclude \"runtime*\" --exclude \"runlog\" --exclude \"*.o[0-9]*\" --exclude \"*.e[0-9]*\" RUN_$runname/ RUN_$descriptor/;
              ln -s RUN_$descriptor RUN1_$descriptor", $return);
//        foreach ($return as $tmp) {
//           echo "$tmp<br>\n";
//        }
//        echo "<br>\n";
        echo "<h3>Run directory created successfully:  Descriptor = $descriptor</h3>\n";
        echo "<h3><a href=\"../run_setup/setup.php?runname=$descriptor\">Setup Run</a> to continue.</h3>\n";
     } else {
        echo "<h2>No code specified, try again.</h2>\n";
     }
   } else {
     echo "<h2>No descriptor specified, try again.</h2>\n";
   }
?>

<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>

</body>
</html>
