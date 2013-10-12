<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<html>
<head>
<title>SWMF GUI</title>
</head>
<body>

<?php include("site_header.php"); ?>
<?php include("site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI</h1>
<br>

<?php
     $return = "";
     Exec("whoami", $return);
     foreach ($return as $tmp) {
        echo "<h2>Running GUI as user: $tmp</h2>";
     }
 ?>

<br><br><br>
<h3>Select an item from the menu on the left.</h3>
<br>

<?php include("site_footer.php"); ?>

</body>
</html>
