<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>
<?php if (isset($_SERVER['PHP_SELF'])) { $scriptname = basename($_SERVER['PHP_SELF']); } ?>

<html>
<head>
<title>SWMF GUI: Create New Run Directory</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Create New Run Directory</h1>

<br>
<h4>Select an existing compiled code to attach to your run.<br>
  Then create a run descriptor (no spaces) to describe the run,
  like Bz_flip.</h4>
<br>

<?php
   $codes = array(); 
   $tmpdir = opendir( "$HOMEDIR/codes" );
   while( $code = readdir( $tmpdir ) ) {
     if (eregi("CODE3_", $code)) {
       $explode = explode("_", $code, 2);
       $name = $explode[1];
       $codes[] = $name; 
     }
   }
   sort($codes);
?>

<FORM TARGET="" ACTION="<?php echo $scriptname; ?>">
<SELECT NAME=codename>
   <OPTION VALUE="">- Compiled SWMF Codes&nbsp;</OPTION>
   <?php
      foreach ($codes as $code) {
         echo "   <OPTION VALUE=\"$code\"";
         if($code == $codename) { echo " SELECTED"; }
         echo ">$code</OPTION>\n";
      }
   ?>
</SELECT>
<INPUT TYPE=submit VALUE="Select"></FORM>
</FORM>
<BR CLEAR=ALL>

<?php if($codename) {
   echo "
<FORM TARGET=\"\" ACTION=\"createrundescriptor.php\">
  <INPUT TYPE=hidden name=codename value=\"$codename\">
  <LABEL>Run Descriptor: </LABEL><INPUT type=text name=\"descriptor\">
  <INPUT TYPE=submit VALUE=\"Create New Run\"> and attach to code \"$codename\"</FORM>
<BR CLEAR=ALL>\n
\n
<h2>Installation Report For \"$codename\"</h2>\n
<pre>\n
./Config.pl -show\n";
   $return = "";
   Exec("cd ../codes/CODE_$codename; ./Config.pl -show", $return);
   foreach ($return as $tmp) {
      echo "$tmp\n";
   }
   echo "\nContent of bin directory:\n";
   $return = "";
   Exec("cd ../codes/CODE_$codename; ls -l bin", $return);
   foreach ($return as $tmp) {
      echo "$tmp\n";
   }
   echo "
</pre>\n
<br>\n
<a href=\"viewlog.php?codename=$codename&logfile=command.log\" TARGET=\"_log\">\n
 View issued command log</a><br>\n
<a href=\"viewlog.php?codename=$codename&logfile=compile.log\" TARGET=\"_log\">\n
 View compile log</a><br>\n
<BR CLEAR=ALL>\n";
   }
?>

<?php include("../site_footer.php"); ?>

</body>
</html>
