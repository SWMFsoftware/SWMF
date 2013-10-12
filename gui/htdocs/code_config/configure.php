<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // If set, use passed value, otherwise use empty string
  $parameter = array('codename');
  foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>

<html>
<head>
<title>SWMF GUI: Configure Code</title>
</head>
<body>

<?php include("../site_header.php"); ?>
<?php include("../site_sidemenu.php"); ?>

<div id="main" width="100%">

<h1>SWMF GUI: Configure Code "<?php echo "$codename" ?>"</h1>
<BR CLEAR=ALL>

<?php
   if(is_file("../codes/CODE1_$codename/Makefile.conf")) { $ISinstalled = "1"; } else { $ISinstalled = "0"; }
   if($ISinstalled) {
	  // Configure code ======================================================
     $return = "";
     Exec("cd ../codes/CODE1_$codename; ./Config.pl -s", $return);
     echo "
<h2>Configuration Summary</h2>
<div class=\"indent\">
<pre>";
     foreach ($return as $tmp) { echo "$tmp\n"; }
     echo "</pre>
<a href=\"viewlog.php?codename=$codename&logfile=command.log\" TARGET=\"_log\">
 View issued command log</a><br>
</div><BR CLEAR=ALL>

<hr>
<h2>Configuration Options</h2>
<br>

<div class=\"indent\">
<FORM TARGET=\"\" ACTION=\"doconfigure.php\">
  <INPUT TYPE=hidden name=codename value=\"$codename\">
  <INPUT TYPE=hidden name=install value=\"-uninstall\">
Uninstall to throw out configuration changes and start from scratch.
  <INPUT TYPE=submit VALUE=\"Uninstall\">
</FORM>

<FORM TARGET=\"\" ACTION=\"doconfigure.php\">
  <INPUT TYPE=hidden name=codename value=\"$codename\">
  <TABLE ALIGN=\"LEFT\" BORDER=0 CELLPADDING=2 CELLSPACING=0>
     ";
     $return = "";
     Exec("cd ../codes/CODE1_$codename; ./Config.pl -l", $return);
     foreach ($return as $tmp) {
       $pos = strpos($tmp, "/");
       if ($pos === false) {
         // nothing, try again
       } else {
         // process for component/model pair
         $chunks = explode("/", $tmp, 2);
         $model = $chunks[0];
         $other = $chunks[1];
         $other2 = "";
         $pos = strpos($other, " ");
         if ($pos === false) {
           // last one
           $choice = $other;
         } else {
           // more, get this one and get ready for next
           $chunks = explode(" ", $other, 2);
           $choice = $chunks[0];
           $other2  = $chunks[1];
         }
         echo "
  <TR>
    <TD>$model:</TD>
    <TD><INPUT TYPE=radio name=$model value=\"\" CHECKED>$choice</TD>\n";
         while($other2) {
           $pos = strpos($other2, "/");
           if ($pos === false) {
             // nothing left
             $other2 = "";
           } else {
             $chunks = explode("/", $other2, 2);      
             $other2 = $chunks[1];
             $pos = strpos($other2, " ");
             if ($pos === false) {
               // last one
               $choice = $other2;
               $other2 = "";
             } else {
               $chunks = explode(" ", $other2, 2);      
               $choice = $chunks[0];
               $other2 = $chunks[1];
             }
             echo "    <TD><INPUT TYPE=radio name=$model value=\"$model/$choice\">$choice</TD>\n";
           }
         }
         echo "  </TR>\n";
       }
     }
     echo "
  </TABLE>
  <br CLEAR=ALL><br>
  <LABEL>Extra configure options: </LABEL><INPUT type=text name=\"options\">
  <a href=\"showhelp.php?codename=$codename\" TARGET=\"_help\">View Config.pl help for options</a><br><br>
  <INPUT TYPE=submit VALUE=\"Configure\">
</FORM>
</div><BR CLEAR=ALL>
     ";
     $return = "";
     Exec("cd ../codes/CODE1_$codename; ./Config.pl -s", $return);
     $grids = array();
     foreach ($return as $tmp) {
       if (eregi("grid:", $tmp)) { $grids[] = $tmp; }
     }
     if($grids) {
       echo "
<hr>
<h2>Set Component Grid</h2>
<br>
<div class=\"indent\">
       ";
       foreach ($grids as $grid) {
         $chunks = explode("/", $grid, 2);
         $model = trim($chunks[0]);
         $grid = $chunks[1];
			echo "
<FORM TARGET=\"\" ACTION=\"doconfigure.php\">
  <INPUT TYPE=hidden name=codename value=\"$codename\">
  <INPUT TYPE=hidden name=grid value=\"$model\">
  <INPUT TYPE=submit VALUE=\"Set grid for $model\">
         ";
         $chunks = explode("grid:", $grid, 2);
         $grid = $chunks[1];
         $i = "0";
			while(ereg(",", $grid)) {
           $chunks = explode(",", $grid, 2);
           $value = trim($chunks[0]);
           $grid = $chunks[1];
           $i += 1;
           echo "<INPUT type=text size=5 name=grid$i value=\"$value\">";
			}
         $value = trim($grid);
         $i += 1;
         echo "<INPUT type=text size=5 name=grid$i value=\"$value\">";
         echo "
  </FORM>
         ";
		 }
       echo "
</div>
<BR CLEAR=ALL>
       ";
     }
     echo "
<hr>
<h2>I'm done configuring the code, prepare for compiling.</h2>
<br>
<div class=\"indent\">
<FORM TARGET=\"\" ACTION=\"movecode.php\">
<INPUT TYPE=hidden name=codename value=\"$codename\">
<INPUT TYPE=submit VALUE=\"DONE Configuring\">
</FORM>
</div>
     ";
   } else {
	  // Install code ======================================================
     echo "
<h2>The code must be installed before it is configured.</h2>
<BR CLEAR=ALL>

<hr>
<h2>Installation Options</h2>
<br>

<FORM TARGET=\"\" ACTION=\"doconfigure.php\">
  <INPUT TYPE=hidden name=codename value=\"$codename\">
     ";
     $return = "";
     Exec("uname", $return);
     $OS = $return[0];
     echo "
  <LABEL>OS Type: $OS </LABEL><br>
  <TABLE ALIGN=\"LEFT\" BORDER=0 CELLPADDING=2 CELLSPACING=0>
  <TR>
    <TD>Compiler:</TD>
	 <TD>
     ";
     $return = "";
     Exec("cd ../codes/CODE1_$codename/share/build; ls Makefile.${OS}*", $return);
     $CHKD = "0";
     foreach ($return as $tmp) {
       $chunks = explode("$OS.", $tmp, 2);
       $makefile = $chunks[1];
       $return2 = "";
       Exec("grep -i COMPILE.f90 ../codes/CODE1_$codename/share/build/$tmp | grep PATH_F", $return2);
       foreach ($return2 as $tmp2) {
         $chunks = explode("PATH_F}", $tmp2, 2);
         $tmp3 = $chunks[1];
	 $pos = strpos($tmp3, " ");
	 if ($pos === false) {
           $compiler = $tmp3;
	 } else {
           $chunks = explode(" ", $tmp3, 2);
           $compiler = $chunks[0];
	 }
       }
       $return2 = "";
       Exec("which $compiler | grep -v 'not found'", $return2);
       if($return2) {
         if (! $CHKD) {
	   $CHKD = "1";
           echo "      <INPUT TYPE=radio name=compiler value=\"$makefile\" CHECKED>$compiler\n";
         } else {
           echo "      <INPUT TYPE=radio name=compiler value=\"$makefile\">$compiler\n";
         }
       }
     }
     echo "
    </TD>
  </TR>
  <br>
  <TR>
    <TD>Named MPI version:</TD>
	 <TD>
     ";
     $return = "";
     Exec("cd ../codes/CODE1_$codename/share/include; ls mpif90_${OS}*", $return);
     //print_r($return);
     $CHKD = "0";
     foreach ($return as $tmp) {
       $chunks = explode(".h", $tmp, 2);
       $tmp = $chunks[0];
       $chunks = explode("${OS}_", $tmp, 2);
       $mpiversion = $chunks[1];
       if (! $CHKD) {
	   $CHKD = "1";
         echo "      <INPUT TYPE=radio name=mpiversion value=\"$mpiversion\" CHECKED>$mpiversion\n";
       } else {
         echo "      <INPUT TYPE=radio name=mpiversion value=\"$mpiversion\">$mpiversion\n";
       }
     }
     echo "
    </TD>
  </TR>
  </TABLE>
  <INPUT TYPE=hidden name=install value=\"-install\">
  <BR CLEAR=ALL>
  <BR CLEAR=ALL>
  <INPUT TYPE=submit VALUE=\"Install\">
</FORM>";
	}
?>

<BR CLEAR=ALL>

<?php include("../site_footer.php"); ?>
