<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // Set posted value
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = $_POST[$name]; } ?>
<?php    // Set posted value
  $parameter = array('name', 'first', 'last', 'stride'); $i = 0;
  while($i < 20) {
    foreach($parameter as $name) {
      $value = "$name$i"; $$value = $_POST[$value]; }
    $i += 1; } ?>

<html>
<head>
<title>SWMF GUI: Save Files</title>

<?php
  $testfile = "../runs/RUN_$runname/LAYOUT.in";
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
  $data = "#COMPONENTMAP\n";
  $i = "1";
  $value = "name$i"; $nam = $$value;
  while($nam) {
    $value = "first$i";  $fir = $$value;
    $value = "last$i";   $las = $$value;
    $value = "stride$i"; $str = $$value;
    $data .= "$nam  $fir  $las $str\n";
    $i += 1;
    $value = "name$i"; $nam = $$value;
  }
  $data .= "#END\n";
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
