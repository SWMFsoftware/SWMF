<?php    // Set posted value
  // Required values
  $parameter = array('plotapplication', 'cmp', 'plottype', 'plotfile', 'filedir');
  foreach($parameter as $name) { $$name = $_POST[$name]; }

  if($plotapplication == "tecplot")    { include("plot_3Dplt.php"); }
  if($plotapplication == "tecplot1D")  { include("plot_1Dlog.php"); }
  if($plotapplication == "tecplot2D")  { include("plot_2Dplt.php"); }
  if($plotapplication == "tecplotLOS") { include("plot_LOS.php"); }
  if($plotapplication == "IEidl")      { include("plot_IEidl.php"); }
  if($plotapplication == "2Didl")      { include("plot_2Didl.php"); }
 ?>

<?php include("plotfunctions.php");?>
<?php include("paths.php");?>

<?php

   echo "
<html>
<head><title>Plot Results: $plotfile</title></head>
<body>
   ";

   // CHECK FOR ERRORS
   $errors = plot0();
   if ($errors) {
     echo "<CENTER>ERROR, incomplete information given.  Fix and try again.</CENTER><br>\n"; exit();
   }

   $plotfileclip = $plotfile;
   $pieces = explode("$plotextension", $plotfile);
   $plotfileclip = $pieces[0];
   $tmpdir = `date +tp%H%M%Stp`;
   $tmpdir = trim($tmpdir);
   $imagedir = "$filedir/../images/${cmp}_$plottype";
   $batchdir = "$filedir/../images/batch";
   if (! is_dir("$batchdir")) { Exec("cd $filedir/../images; mkdir batch"); }

   // LOAD 1ST CUSTOM CODE BLOCK TO CREATE SCRIPTS
   plot1();

   $reusemacro = "0";
   $files = array(); 
   $dir2 = opendir($imagedir);
   while( $file = readdir( $dir2 ) ) {
     if (eregi("batch", $file) && eregi("$macroextension", $file)) {
       $files[] = $file; 
       $return = "";
       Exec("diff $imagedir/$file $batchdir/$tmpdir/batch$macroextension;", $return);
       $diffs = count($return);
       if ($diffs == 0) {
         $pieces = explode(".", $file);
         $file = $pieces[0];
         $pieces = explode("-", $file);
         $reusenumber = $pieces[1];
         $reusemacro = "1";
       }
     }
   }
   sort($files);
   $number = count($files);
   $file = $files[$number-1];
   $pieces = explode(".", $file);
   $file = $pieces[0];
   $pieces = explode("-", $file);
   $number = $pieces[1];
   $number += 1;
   if($number < 10) { $number = "00$number"; } elseif($number < 100) { $number = "0$number"; }
   if ($reusemacro) {
     $number = $reusenumber;
   }
   $file1 = "${plotfileclip}-${number}.ps";
   $file2 = "${plotfileclip}-${number}.png";
   $file3 = "batch-${number}$macroextension";

   $fileexists = "0";
   if ($reusemacro) {
     $dir2 = opendir($imagedir);
     while( $file = readdir( $dir2 ) ) {
       if (eregi("$plotfileclip-$number", $file)) {
         $fileexists = "1";
       }
     }
   } else {
     Exec("cd $batchdir/$tmpdir;
           cp batch$macroextension ../../${cmp}_$plottype/$file3");
     // LOAD 2ND CUSTOM CODE BLOCK TO SAVE DEFAULTS FILE
     plot2();
   }

   $decodedplotfile = decodeFilename($plotfile);
   echo "
<CENTER>File=$plotfile, $decodedplotfile, Style=$number</CENTER><br>
   ";

   if (!($fileexists)) {
     echo "
<CENTER>Please wait while image is rendered ...</CENTER>
     ";
     ob_flush();flush();
     $time1 = time();
     set_time_limit(0);

     // LOAD 3RD CUSTOM CODE BLOCK TO RUN SCRIPT
     plot3();

     Exec("cd $batchdir/$tmpdir;
           $gs -sDEVICE=png16m -sOutputFile=tmp.png -dNOPAUSE -q -dBATCH $file1;
           $convert +antialias -trim tmp.png $file2;
           cp $file1 ../../${cmp}_$plottype/;
           cp $file2 ../../${cmp}_$plottype/");
     $time2 = time();
     $timedif = $time2-$time1;
     echo "
<CENTER>processing time: $timedif seconds</CENTER><br>
     ";
   } else {
     echo "
<CENTER>Using previously created plot.</CENTER><br>
     ";
   }
   Exec("cd $batchdir;
         rm -rf $tmpdir");

   echo "
<P><CENTER>
<A HREF=\"$filedir/../images/${cmp}_$plottype/$file1\"><IMG SRC=\"$filedir/../images/${cmp}_$plottype/$file2\" BORDER=0></A>
</CENTER></P>
<BR CLEAR=ALL>
   ";

   // LOAD 4TH CUSTOM CODE BLOCK TO DISPLAY DOWNLOAD FOR SELF-PLOT INSTRUCTIONS
   plot4();

?>
