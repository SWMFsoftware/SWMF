<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php
  include("paths.php");

  $plotfileclip = $plotfile;
  $pieces = explode("$plotextension", $plotfile);
  $plotfileclip = $pieces[0];
  $tmpdir = `date +tp%H%M%Stp`;
  $tmpdir = trim($tmpdir);
  $imagedir = "$filedir/../images/${cmp}_$plottype";
  $batchdir = "$filedir/../images/batch";
  if (! is_dir("$batchdir")) { Exec("cd $filedir/../images; mkdir batch"); }

  $file1 = "${plotfileclip}-${number}.ps";
  $file2 = "${plotfileclip}-${number}.png";
  $file3 = "batch-${number}$macroextension";

  $fileexists = "0";
  $dirTMP = opendir($imagedir);
  while( $file = readdir( $dirTMP ) ) {
    if (eregi("$plotfileclip-$number", $file)) {
      $fileexists = "1";
    }
  }

  if (!($fileexists)) {
    $decodedplotfile = decodeFilename($plotfile);

    Exec("mkdir $batchdir/$tmpdir");

    include("quick_${loadfile}.php");

    Exec("cd $batchdir/$tmpdir;
          echo '#!/bin/sh' > runscript.sh;
          echo '' >> runscript.sh;
          echo '${gs} -sDEVICE=png16m -sOutputFile=tmp.png -dNOPAUSE -q -dBATCH ${file1}' >> runscript.sh;
          echo '${convert} +antialias -trim tmp.png ${file2}' >> runscript.sh;
          echo 'cp ${file1} ../../${cmp}_${plottype}/' >> runscript.sh;
          echo 'cp ${file2} ../../${cmp}_${plottype}/' >> runscript.sh;
          echo '' >> runscript.sh;
          chmod 755 runscript.sh");

    Exec("cd $batchdir/$tmpdir; ./batchscript.sh");
    Exec("cd $batchdir/$tmpdir; ./runscript.sh");

    Exec("rm -rf $batchdir/$tmpdir");
  }
?>
