<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<html>
<head>
<title>Upload File</title>
<META HTTP-EQUIV=REFRESH CONTENT="2;URL=index.php">
</head>
<body>

<h3>This log window will clear shortly, please wait ...</h3>

<?php
   set_time_limit(0);
   $target_path = "";
   $target_path = $target_path . basename($_FILES['uploadedfile']['name']);
   if(move_uploaded_file($_FILES['uploadedfile']['tmp_name'], $target_path)) {
      echo "The file " . basename($_FILES['uploadefile']['name']) . " has been uploaded";
   } else {
      echo "There was an error uploading the file, please try again!";
   }
?>

</body>
