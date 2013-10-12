<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('plotname');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Manage Plots</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Manage Plots</h1>
<BR CLEAR=ALL>

<?php
  if($plotname) {
    echo "
<br>
<h2><a href=\"/plot_manage/index.php\">Return</a> to manage plots without changes.</h2>
<br><br><br>
<h2>Plot selected: $plotname</h2>
<br>
<h2><a href=\"/plot_manage/remove.php?plotname=$plotname\">REMOVE</a> (THERE IS NO UNDO ONCE CLICKED!)</h2>
<br><br>
<h2><a href=\"/plot_manage/clean.php?plotname=$plotname\">CLEAN</a> Delete all created plots and styles only.</h2>
    ";

  } else {
    echo "
<h2>View plot directories for deletion:</h2><br>
<div class=\"indent\">
<table>
    ";
    $tmpdir = opendir( "$HOMEDIR/plots" );
    while( $plot = readdir( $tmpdir ) ) {
      if (ereg("PLOT_", $plot)) { $plots[] = $plot; }
    }
    if($plots){ sort($plots);
      foreach ($plots as $plot) {
        $explode = explode("_", $plot, 2); $name = trim($explode[1]);
        echo "
<tr>
<td><b><a href=\"/plot_manage/index.php?plotname=$name\">$name</a></b></td>
<td>&nbsp&nbsp</td>
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
