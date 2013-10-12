<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // Set posted value
  $parameter = array('tdir');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; }
  if($tdir) {
    $time1 = time();
    set_time_limit(0);
    Exec("cd $tdir; ./batchscript.sh");
    Exec("cd $tdir; ./runscript.sh");
    $time2 = time();
    $timedif = $time2-$time1;
    Exec("cd $tdir;
          echo '<br><hr><CENTER>processing time: ${timedif} seconds</CENTER><br>' >> page.html");

    $fname = "$tdir/page.html";
    $fh = fopen($fname, "rt");
    $fcontent = fread($fh, filesize($fname));
    fclose($fh); 
    echo $fcontent;

    Exec("rm -rf $tdir");
    exit();
  }

  // Required values
  $parameter = array('cmp', 'plottype', 'plotfile', 'filedir');
  foreach($parameter as $name) { $$name = ""; $$name = $_POST[$name]; }

  if($cmp == "") {
     echo "
<html>
<head><title>ERROR!</title></head>
<body><CENTER>ERROR, no component given.  This page can't be bookmarked.</CENTER><br></body>
</html>";
     exit();
  }

  // Load variables and functions
  include("${filedir}/../images/${cmp}_${plottype}/defaultsBASE.php");
  include("plot_${loadfile}.php");
  include("plotfunctions.php");
  include("paths.php");

  // CHECK FOR PLOT SPECIFIC ERRORS
  $errors = plot0();
  if ($errors) {
    echo "
<html>
<head><title>ERROR!</title></head>
<body><CENTER>ERROR, incomplete information given.  Fix and try again.</CENTER><br></body>
</html>";
    exit();
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

  if (!($fileexists)) {
    Exec("cd $batchdir/$tmpdir;
          echo '<html>' >> page.html;
          echo '<head>' >> page.html;
          echo '<title>Plot Results: ${plotfile}</title>' >> page.html;
          echo '</head>' >> page.html;
          echo '<body>' >> page.html;
          echo '<CENTER>File=${plotfile}, ${decodedplotfile}, Style=${number}</CENTER><br>' >> page.html;
          echo '' >> page.html");

    // LOAD 3RD CUSTOM CODE BLOCK TO RUN SCRIPT
    plot3();

    Exec("cd $batchdir/$tmpdir;
          echo '#!/bin/sh' > runscript.sh;
          echo '' >> runscript.sh;
          echo '${gs} -sDEVICE=png16m -sOutputFile=tmp.png -dNOPAUSE -q -dBATCH ${file1}' >> runscript.sh;
          echo '${convert} +antialias -trim tmp.png ${file2}' >> runscript.sh;
          echo 'cp ${file1} ../../${cmp}_${plottype}/' >> runscript.sh;
          echo 'cp ${file2} ../../${cmp}_${plottype}/' >> runscript.sh;
          echo '' >> runscript.sh;
          chmod 755 runscript.sh");
  } else {
    Exec("cd $batchdir/$tmpdir;
          echo '<html>' >> page.html;
          echo '<head>' >> page.html;
          echo '<title>Plot Results: ${plotfile}</title>' >> page.html;
          echo '</head>' >> page.html;
          echo '<body>' >> page.html;
          echo '<CENTER>File=${plotfile}, ${decodedplotfile}, Style=${number}</CENTER><br>' >> page.html;
          echo '' >> page.html");
  }

  Exec("cd $batchdir/$tmpdir;
        echo '<P><CENTER>' >> page.html;
        echo '<A HREF=\"${filedir}/../images/${cmp}_${plottype}/${file1}\"><IMG SRC=\"${filedir}/../images/${cmp}_${plottype}/${file2}\" BORDER=0></A>' >> page.html;
        echo '</CENTER></P>' >> page.html;
        echo '<BR CLEAR=ALL>' >> page.html;
        echo '' >> page.html");

  // LOAD 4TH CUSTOM CODE BLOCK TO DISPLAY DOWNLOAD FOR SELF-PLOT INSTRUCTIONS
  plot4();

  if (!($fileexists)) {
    echo "
<html>
<head>
<title>Plot Results: ... rendering ...</title>
<META HTTP-EQUIV=\"refresh\" content=\"0;URL=createplot.php?tdir=${batchdir}/${tmpdir}\">
</head>
<body>
<br><center><img src=\"../images/pleasewait.gif\"></center>
</body>
</html>
    ";
  } else {
    Exec("cd $batchdir/$tmpdir;
          echo '<br><hr><CENTER>Using previously created plot.</CENTER><br>' >> page.html;
          echo '' >> page.html");

    $fname = "$batchdir/$tmpdir/page.html";
    $fh = fopen($fname, "rt");
    $fcontent = fread($fh, filesize($fname));
    fclose($fh); 
    echo $fcontent;

    Exec("rm -rf $batchdir/$tmpdir");
  }

?>
