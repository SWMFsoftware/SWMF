<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<html>
<head>
<title>SWMF GUI: Slice-N-Dice For Perry</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Slice-N-Dice For Perry</h1>
<BR CLEAR=ALL>

<h2>3D files and status: (NOTE: Now png file output)</h2><br>
<div class=\"indent\">
<table>

<?php
    $tmpdir = opendir( "." );
    while( $file = readdir( $tmpdir ) ) {
      if (ereg(".plt", $file)) { $files[] = $file; }
    }
    if($files){
      sort($files);
      foreach ($files as $file) {
        $status = "";
        $explode = explode(".plt", $file, 2); $name = $explode[0];
        if(is_file("${name}.tgz")) {
          $status = "processing completed. <a href=\"${name}.tgz\">download tarball</a>";
        } elseif(is_file("$name/README")) {
          $status = "processing ...";
        } else {
          $status = "<a href=\"process.php?file=$name\">process</a>";
        }
        echo "
<tr>
<td><b>$file</b></td>
<td>&nbsp&nbsp</td><td>$status</td>
</tr>
        ";
      }
    }
    echo "
</table>
</div>
    ";
?>
<BR CLEAR=ALL>

<h2>Upload new 3D file:</h2><br>
<FORM TARGET="" METHOD="post" ACTION="upload.php" ENCTYPE="multipart/form-data">
  <INPUT TYPE=hidden name=MAX_FILE_SIZE VALUE=500000000>
  <INPUT TYPE=file NAME=uploadedfile><br>
  <INPUT TYPE=submit VALUE="Upload"> a new file to run directory (500mb limit) &nbsp&nbsp&nbsp
</FORM>

<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>

</body>
</html>
