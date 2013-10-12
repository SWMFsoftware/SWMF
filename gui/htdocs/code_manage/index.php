<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Manage Codes</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Manage Codes</h1>
<BR CLEAR=ALL>

<?php
  if($codename) {
    echo "
<br>
<h2><a href=\"/code_manage/index.php\">Return</a> to manage codes without changes.</h2>
<br><br><br>
<h2>Code selected: $codename</h2>
<br>
<h2><a href=\"/code_manage/remove.php?codename=$codename\">REMOVE</a> (THERE IS NO UNDO ONCE CLICKED!)</h2>
    ";

  } else {
    echo "
<h2>View codes for deletion:</h2><br>
<div class=\"indent\">
<table>
    ";
    $codes = "";
    $tmpdir = opendir( "$HOMEDIR/codes" );
    while( $code = readdir( $tmpdir ) ) {
      if (ereg("CODE_", $code)) { $codes[] = $code; }
    }
    if($codes){ sort($codes);
      foreach ($codes as $code) {
        $explode = explode("_", $code, 2); $name = trim($explode[1]);
        $inst = "Installed"; $conf = ""; $comp = "";
        if (is_file("$HOMEDIR/codes/CODE2_${name}/Makefile")) { $conf = "/ Configured"; }
        if (is_file("$HOMEDIR/codes/CODE3_${name}/Makefile")) { $conf = "/ Configured"; $comp = "/ Compiled"; }
        echo "
<tr>
<td><b><a href=\"/code_manage/index.php?codename=$name\">$name</a></b></td>
<td>&nbsp&nbsp</td><td>$inst</td><td>$conf</td><td>$comp</td>
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
