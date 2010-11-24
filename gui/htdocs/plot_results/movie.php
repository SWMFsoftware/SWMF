<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('runname', 'cmp', 'plottype', 'plotstyle', 'wait');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; }

  if (($wait)) {
    echo "
<html>
<head>
<title>View Movie: ... rendering ...</title>
<META HTTP-EQUIV=\"refresh\" content=\"0;URL=movie.php?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=$plotstyle\">
</head>
<body>
<br><center><img src=\"../images/pleasewait.gif\"></center>
</body>
</html>
    ";
    exit();
  }

  $time1 = time();
  set_time_limit(0);

  include("plotfunctions.php");

  $runpath = "../plots/PLOT_$runname";
  $filedir = "$runpath/$cmp";

  include("$runpath/images/${cmp}_$plottype/defaultsBASE.php");
  $number = "$plotstyle";

  $plotfile = "";
  $choice = "";
  include("$runpath/images/${cmp}_$plottype/defaults-${number}.php");
  if ("$choice" == "") {
    $plotfilelist = GetPlotList("$runpath/$cmp");
  } else {
    $plotfilelist = GetXYZPlotList("$runpath/$cmp", "$choice");
  }

  $newfiles = "0";
  $countfiles = count($plotfilelist);
  if($countfiles > 0) {
    foreach ($plotfilelist as $plotfile) {
      include("makeplot.php");
      if (!($fileexists)) { $newfiles = "1"; }
    }
  }
  if ($newfiles || !(is_file("$runpath/images/${cmp}_${plottype}/animation-${number}.mp4"))) {
    Exec("cd $runpath/images/${cmp}_${plottype};
          ../../../../plot_results/movie.sh $number");
  }

  $time2 = time();
  $timedif = $time2-$time1;
?>
<html>
<head>
<title>View Movie</title>
</head>
<body>
<?php echo "
<h2>Movie of files for run $runname - $cmp - $plottype - $plotstyle </h2>
<h3>$countfiles frames in animation, proccessing time ${timedif} seconds</h3>
<center>
<EMBED SRC=\"${filedir}/../images/${cmp}_${plottype}/animation-${number}.mp4\"
 WIDTH=\"640\" HEIGHT=\"480\"
 AUTOPLAY=\"true\" CONTROLLER=\"true\"
 PLUGINSPAGE=\"http://www.apple.com/quicktime/download/\">
</EMBED>
</center>
"; ?>

</body>
</html>
