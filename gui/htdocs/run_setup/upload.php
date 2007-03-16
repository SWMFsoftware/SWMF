<?php    // Set posted value
  $parameter = array('runname');
  foreach($parameter as $name) { $$name = $_POST[$name]; } ?>

<html>
<head>
<title>SWMF GUI: Upload File</title>
<META HTTP-EQUIV=REFRESH CONTENT="10;URL=setup.php?runname=<?php echo $runname ?>">
</head>
<body>

<h3>This log window will clear shortly, please wait ...</h3>

<?php
   $target_path = "../runs/RUN_$runname/";
   $target_path = $target_path . basename($_FILES['uploadedfile']['name']);
   if(move_uploaded_file($_FILES['uploadedfile']['tmp_name'], $target_path)) {
      echo "The file " . basename($_FILES['uploadefile']['name']) . " has been uploaded";
   } else {
      echo "There was an error uploading the file, please try again!";
   }
?>

</body>
