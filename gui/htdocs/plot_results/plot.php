<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Visualize Run Output</title>
</head>
<body>

<?php
 include("../site_header.php");
 include("../site_sidemenu.php");

 include("plotheader.php");
 include("plotfunctions.php");
?>

<div id="righticonbar">
<br><br>
<?php
 $runpath = "../plots/PLOT_$runname";
 $margin = 150;
 include("plotright.php");
?>
</div>

<div id="main2" width="100%">

<?php
 include("plotcenter.php");
 include("../site_footer.php");
?>

</body>
</html>
