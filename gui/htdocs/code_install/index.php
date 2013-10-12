<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<html>
<head>
<title>SWMF GUI: Create New Code Directory</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Create New Code Directory</h1>

<br>
<h4>Create a code descriptor (no spaces) to describe your
  configuration, like Darren_GM-IE-IM_Earth.</h4>
<br>
<TABLE ALIGN="LEFT" BORDER=0 CELLPADDING=0 CELLSPACING=0><TR>
  <TD><FORM TARGET="" METHOD="post" ACTION="createcodedescriptor.php">
      <LABEL>Code Descriptor: </LABEL><INPUT type=text name="descriptor">
      <INPUT TYPE=submit VALUE="Submit"></FORM></TD>
</TR></TABLE>
<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>

</body>
</html>
