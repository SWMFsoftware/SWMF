<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Compile Code</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Compile Code "<?php echo "$codename" ?>"</h1>
<BR CLEAR=ALL>

<h2>Current Code Summary</h2>
<pre>
<?php
   echo "./Config.pl -s\n";
   $return = "";
   Exec("cd ../codes/CODE2_$codename; ./Config.pl -s", $return);
   foreach ($return as $tmp) {
      echo "$tmp\n";
   }
   echo "\nContent of bin directory:\n";
   $return = "";
   Exec("cd ../codes/CODE2_$codename; ls -l bin", $return);
   foreach ($return as $tmp) {
      echo "$tmp\n";
   }
?>
</pre>
<br>
View <a href="viewlog.php?codename=<?php echo $codename ?>&logfile=command.log" TARGET="_log">
 issued command log</a><br>
<BR CLEAR=ALL>

<hr>
<h2>Compile Code</h2>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp This may take a few minutes!</h3>
<br>
<FORM TARGET="" ACTION="docompile.php">
  <INPUT TYPE=hidden name=codename value="<?php echo $codename ?>">
  <INPUT TYPE=hidden name=options value="ALL">
  <INPUT TYPE=submit VALUE="make ALL">
  <LABEL>(most common compile option)</LABEL>
</FORM>
<FORM TARGET="" ACTION="docompile.php">
  <INPUT TYPE=hidden name=codename value="<?php echo $codename ?>">
  <INPUT TYPE=submit VALUE="make">
  <LABEL>with options: </LABEL><INPUT type=text name="options">
  <a href="showmake.php?codename=<?php echo $codename ?>" TARGET="_make">View make options</a><br>
</FORM>
<BR CLEAR=ALL>
<iframe width=105% height=230 frameborder=1 src="taillog.php?codename=<?php echo $codename; ?>"></iframe>
<br><br>

<hr>
<?php
  $return = "";
  Exec("ps uxc | grep make | grep -v grep", $return);
  if($return) {
    echo "
<h2>Compile in progress ...</h2>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp When compilation has finished, refresh page for further options.</h3>
<br>
    ";
  } else {
    echo "
<h2>Make Code Available For Runs</h2>
<h3>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp You cannot go back and compile after moving code.</h3>
<br>
<FORM TARGET=\"\" ACTION=\"movecode.php\">
  <INPUT TYPE=hidden name=codename value=\"$codename\">
  <INPUT TYPE=submit VALUE=\"Move Code\"></FORM>
    ";
  }
?>

<?php include("../site_footer.php"); ?>

</body>
</html>
