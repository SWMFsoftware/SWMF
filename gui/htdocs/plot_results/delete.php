<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname', 'cmp', 'plottype', 'plotstyle', 'plotfile', 'delete');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<?php
if($delete) {
  echo "
<html><head>
<META HTTP-EQUIV=\"refresh\" content=\"0;URL=plot.php?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=001&plotfile=$plotfile\">
</head><body>
  ";
  Exec("cd ../plots/PLOT_$runname/images/${cmp}_$plottype; rm -f *-${plotstyle}.*");
  echo "</body>";
  exit();
}
?>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<html>
<head>
<title>SWMF GUI: Delete individual plot styles</title>
</head>
<body>

<h1>Delete individual plot styles</h1>
<BR CLEAR=ALL>

<?php
  echo "
<br>
<h2><a href=\"plot.php?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=$plotstyle&plotfile=$plotfile\">Return</a> without changes.</h2>
<br><br>
<h2><a href=\"delete.php?delete=1&runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=$plotstyle&plotfile=$plotfile\">REMOVE ALL</a> (THERE IS NO UNDO ONCE CLICKED!)</h2>
<br><br>
<h2>Style:</h2>
<br>
<div class=\"indent\">
run name = $runname<br>
component = $cmp<br>
plot type = $plottype<br>
plot style = $plotstyle<br>
</div>
<br>
<h2>Existing plots with this style:</h2>
<br>
<div class=\"indent\">
  ";
  $tmpdir = opendir( "../plots/PLOT_$runname/images/${cmp}_$plottype" );
  while( $type = readdir( $tmpdir ) ) {
    if (ereg("-${plotstyle}.png", $type)) {
      echo "$type<br>";
    }
  }
  echo "</div>";
?>
<br><br>

</body>

<?php include("../site_footer.php"); ?>

</body>
</html>
