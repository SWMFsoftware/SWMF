<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('refresh');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; }
?>

<html>
<head>
<title>TEST</title>
<?php
  if(!$refresh) {
    echo "<META HTTP-EQUIV=\"refresh\" content=\"1;URL=pleasewait.php?refresh=1\">";
  }
?>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<?php
  if(!$refresh) {
    echo "<br><center><img src=\"../images/pleasewait.gif\"></center>";
    exit();
  }
?>

<h1>TEST</h1>

<br>
<h4>bla bla bla
  </h4>
<br>

<?php Exec("sleep 15"); ?>

<h4>ble ble ble ble</h4>

<?php include("../site_footer.php"); ?>

</body>
</html>
