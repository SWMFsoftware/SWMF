<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // Set posted value
  $parameter = array('runname', 'filename', 'data');
  foreach($parameter as $name) { $$name = $_POST[$name]; } ?>

<html>
<head>
<title>SWMF GUI: Save Files</title>

<?php
   $testfile = "../runs/RUN_$runname/$filename";
   if(!$fh = fopen($testfile, 'w')) {
     echo "
<META HTTP-EQUIV=REFRESH CONTENT=\"10;URL=setup.php?runname=$runname\">
</head>
<body>
The file $testfile cannot be opened
</body>
</html>
     ";
     exit;
   }
   if(fwrite($fh, $data) === FALSE) {
     echo "
<META HTTP-EQUIV=REFRESH CONTENT=\"10;URL=setup.php?runname=$runname\">
</head>
<body>
The file $testfile cannot be written to
</body>
</html>
     ";
     exit;
   }
   fclose($fh);
   Exec("../../bin/dos2unix.pl $testfile", $return);
?>

<META HTTP-EQUIV=REFRESH CONTENT="0;URL=setup.php?runname=<?php echo $runname ?>">
</head>
<body>
</body>
</html>
