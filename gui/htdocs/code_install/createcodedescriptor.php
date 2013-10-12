<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php $descriptor = $_POST["descriptor"]; ?>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<?php
  // If no descriptor, print warning and exit.
  if (! $descriptor) {
    echo "
<html>
<head>
<title>SWMF GUI: Create New Code Directory</title>
<META HTTP-EQUIV=REFRESH CONTENT=\"10;URL=.\">
</head>
<body>
<h1>SWMF GUI: Create New Code Directory</h1>
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

  // If code already exists, print warning and exit.
  if (is_dir("../codes/CODE_$descriptor")) {
    echo "
<html>
<head>
<title>SWMF GUI: Create New Code Directory</title>
<META HTTP-EQUIV=REFRESH CONTENT=\"10;URL=.\">
</head>
<body>
<h1>SWMF GUI: Create New Code Directory</h1>
<BR CLEAR=ALL>
<h2>Descriptor already exists, try again.</h2>
<BR CLEAR=ALL>
    ";
    include("../site_footer.php");
    echo "
</body></html>
    ";
    exit();
  }

  // Go ahead and create code
  if (! is_dir("../codes")) { Exec("mkdir ../codes"); }
  $return = "";
  Exec("cd ../codes;
        mkdir CODE_$descriptor;
        rsync -avq --exclude gui --exclude 'run*' ../../../. CODE_$descriptor/;
        ln -s CODE_$descriptor CODE1_$descriptor;
	cd CODE_$descriptor;
	./Config.pl -install >& .tmp;
	./Config.pl -uninstall >& .tmp;
	rm -f .tmp", $return);
//  foreach ($return as $tmp) {
//    echo "$tmp<br>\n";
//  }
//  echo "<br>\n";
  echo "
<html>
<head>
<META HTTP-EQUIV=REFRESH CONTENT=\"0;URL=../code_config/configure.php?codename=$descriptor\">
</head>
<body></body>
</html>
  ";
?>
