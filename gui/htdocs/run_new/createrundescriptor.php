<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename', 'descriptor');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<?php
  // If no codename, print warning and exit.
  if (! $codename) {
    echo "
<html>
<head>
<title>SWMF GUI: Create New Run Directory</title>
<META HTTP-EQUIV=REFRESH CONTENT=\"10;URL=.\">
</head>
<body>
<h1>SWMF GUI: Create New Run Directory</h1>
<BR CLEAR=ALL>
<h2>No codename specified, try again.</h2>
<BR CLEAR=ALL>
    ";
    include("../site_footer.php");
    echo "
</body></html>
    ";
    exit();
  }

  // If no descriptor, print warning and exit.
  if (! $descriptor) {
    echo "
<html>
<head>
<title>SWMF GUI: Create New Run Directory</title>
<META HTTP-EQUIV=REFRESH CONTENT=\"10;URL=.?codename=$codename\">
</head>
<body>
<h1>SWMF GUI: Create New Run Directory</h1>
<BR CLEAR=ALL>
<h2>No descriptor specified, try again.</h2>
<BR CLEAR=ALL>
    ";
    include("../site_footer.php");
    echo "
</body></html>
    ";
    exit();
  }

  // If directory exists, print warning and exit.
  if (is_dir("../runs/RUN_$descriptor") || is_dir("../plots/PLOT_$descriptor")) {
    echo "
<html>
<head>
<title>SWMF GUI: Create New Run Directory</title>
<META HTTP-EQUIV=REFRESH CONTENT=\"10;URL=.?codename=$codename\">
</head>
<body>
<h1>SWMF GUI: Create New Run Directory</h1>
<BR CLEAR=ALL>
<h2>Directory already exists, try again.</h2>
<BR CLEAR=ALL>
    ";
    include("../site_footer.php");
    echo "
</body></html>
    ";
    exit();
  }

  if (! is_dir("../runs")) { Exec("mkdir ../runs"); }
  $return = "";
  Exec("cd ../codes/CODE_$codename;
        mkdir RUN_$descriptor;
        rsync -av runCLEAN/ RUN_$descriptor/;
        echo $codename > RUN_$descriptor/.codename;
        mv RUN_$descriptor ../../runs/;
        cd ../../runs;
        ln -s RUN_$descriptor RUN1_$descriptor", $return);
//  foreach ($return as $tmp) {
//    echo "$tmp<br>\n";
//  }
//  echo "<br>\n";
  echo "
<html>
<head>
<META HTTP-EQUIV=REFRESH CONTENT=\"0;URL=../run_setup/setup.php?runname=$descriptor\">
</head>
<body></body>
</html>
  ";
?>
