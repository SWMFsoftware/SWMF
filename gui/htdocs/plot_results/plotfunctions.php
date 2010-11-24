<?php
////////////////////////////////////////////////////////////////////////////////

if (isset($_SERVER['PHP_SELF'])) {  // Name of this file
   $scriptname = basename($_SERVER['PHP_SELF']);
} else {
   $scriptname = 'plot.php';
}

////////////////////////////////////////////////////////////////////////////////


/// Functions ///

function GetFileList($dirname='.') {   // Finds all the figures
// Defaults to all files in alphabetical order.
   $files = array();
   $dir = opendir($dirname);
   while( $file = readdir( $dir ) ) {
     if (eregi("\.jpe?g$", $file) || eregi("\.gif$", $file) || eregi("\.png$", $file)) {
       $files[] = $file; 
     }
   }
   sort($files);
   return $files;
}

function GetFileListByStyle($dirname='.', $string) {   // Finds all the figures of a given style
// Defaults to all files in alphabetical order.
   $files = array();
   $dir = opendir($dirname);
   while( $file = readdir( $dir ) ) {
     if (eregi("\.jpe?g$", $file) || eregi("\.gif$", $file) || eregi("\.png$", $file) || eregi("\.mp4$", $file)) {
       if (eregi("${string}", $file)) {
         $files[] = $file; 
       }
     }
   }
   sort($files);
   return $files;
}

function GetStyleList($dirname='.') {   // Finds all the batch macro files
// Defaults to all files in alphabetical order.
   $files = array();
   $dir = opendir($dirname);
   while( $file = readdir( $dir ) ) {
     if (eregi("batch", $file)) {
       $files[] = $file; 
     }
   }
   sort($files);
   return $files;
} 

function GetPlotList($dirname='.') {   // Finds all the plt files
// Defaults to all files in alphabetical order.
   global $cmp, $plottype;
   $files = array();
   $dir = opendir($dirname);
   while( $file = readdir( $dir ) ) {
     if($cmp == "GM" || $cmp == "SC" || $cmp == "IH") {
       if($plottype == "3Dplt") {
         if (eregi("3d", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "2Dplt") {
         if (eregi("=", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "2Didl") {
         if (eregi("=0", $file) && eregi("\.out$", $file)) { $files[] = $file; } }
       if($plottype == "LOS") {
         if (eregi("los", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "1Dlog") {
         if (eregi("log", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "1Dsat") {
         if (eregi("sat", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
     } elseif($cmp == "IM") {
       if (eregi("2d", $file) && eregi("\.plt$", $file)) { $files[] = $file; }
     } elseif($cmp == "IE") {
       if($plottype == "idl") {
         if (eregi("\.idl$", $file)) { $files[] = $file; } }
     }
   }
   sort($files);
   return $files;
} 

function GetXYZPlotList($dirname='.', $string) {   // Finds all the plt files
// Defaults to all files in alphabetical order.
   global $cmp, $plottype;
   $files = array();
   $dir = opendir($dirname);
   while( $file = readdir( $dir ) ) {
     if($cmp == "GM" || $cmp == "SC" || $cmp == "IH") {
       if($plottype == "3Dplt") {
         if (eregi("3d", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "2Dplt") {
         if (eregi("${string}=", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "2Didl") {
         if (eregi("${string}=", $file) && eregi("\.out$", $file)) { $files[] = $file; } }
       if($plottype == "LOS") {
         if (eregi("_${string}_", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "1Dlog") {
         if (eregi("log", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
       if($plottype == "1Dsat") {
         if (eregi("sat", $file) && eregi("\.plt$", $file)) { $files[] = $file; } }
     } elseif($cmp == "IM") {
       if (eregi("2d", $file) && eregi("\.plt$", $file)) { $files[] = $file; }
     } elseif($cmp == "IE") {
       if($plottype == "idl") {
         if (eregi("\.idl$", $file)) { $files[] = $file; } }
     }
   }
   sort($files);
   return $files;
} 

function GetPlotVariables($filename='.') {   // Extract plot variables from plot file
  global $cmp, $plottype, $runname, $plotfile;
  include("paths.php");
  $vars = array();
  $return = "";
  if($cmp == "IE") {
    Exec("head -30 $filename", $return);
    $found = 0;
    $i = 0;
    foreach ($return as $tmp) {
      if($found) {
        $i++;
        if(eregi("$i", $tmp)) { 
          $chunks = split("$i", $tmp, 2);
          $var = $chunks[1];
          $vars[] = $var;
        } else { $found = 0; }
      }
      if(eregi("VARIABLE LIST", $tmp)) { $found = 1; }
    }
    return $vars;
  }
  if($plottype == "UNUSED_2Didl") { // UNUSED OLD METHOD
    $checkidlexists = "0";
    $dir2 = opendir("../plots/PLOT_$runname");
    while( $file = readdir( $dir2 ) ) {
      if (eregi("checkidl.sh", $file)) { $checkidlexists = "1"; }
    }
    if(! $checkidlexists) {
      Exec("rsync -a BASE/checkidl.sh ../plots/PLOT_$runname/");
    }
    Exec("cd ../plots/PLOT_$runname; ./checkidl.sh $cmp/$plotfile | grep plotvars", $return);
    foreach ($return as $tmp) {
      if(eregi("plotvars", $tmp)) {
        $chunks = split("=", $tmp, 2);
        $tmp = trim($chunks[1]);
        while(eregi(" ", $tmp)) {
          $chunks = split(" ", $tmp, 2);
          $var = $chunks[0];
          $tmp = $chunks[1];
          $vars[] = $var;
        }
        $vars[] = $tmp;
        return $vars;
      }
    }
  }
  if($plottype == "2Didl") {
    $fixendianexists = "0";
    $dir2 = opendir("../plots/PLOT_$runname");
    while( $file = readdir( $dir2 ) ) {
      if (eregi("FixEndian.pl", $file)) { $fixendianexists = "1"; }
    }
    if(! $fixendianexists) {
      Exec("rsync -a BASE/FixEndian.pl ../plots/PLOT_$runname/");
    }
    Exec("cd ../plots/PLOT_$runname; ./FixEndian.pl -i -e=m $cmp/$plotfile", $return);
    // Find how many elements of each type (ndim, neqpar, nw)
    foreach ($return as $tmp) {
      if(eregi("ndim", $tmp)) {
        $chunks = split("ndim=-", $tmp, 2);
        $tmp = trim($chunks[1]);
        $chunks = split(" ", $tmp, 2);
        $ndim = trim($chunks[0]);
        $tmp = trim($chunks[1]);
        $chunks = split("neqpar=", $tmp, 2);
        $tmp = trim($chunks[1]);
        $chunks = split(" ", $tmp, 2);
        $neqpar = trim($chunks[0]);
        $tmp = trim($chunks[1]);
        $chunks = split("nw=", $tmp, 2);
        $nw = trim($chunks[1]);
      }
    }
    // Get var list
    foreach ($return as $tmp) {
      if(eregi("names", $tmp)) {
        $chunks = split("=", $tmp, 2);
        $tmp = trim($chunks[1]);
        $i = 0;
        while(eregi(" ", $tmp)) {
          $i++;
          $chunks = split(" ", $tmp, 2);
          $var = $chunks[0];
          $tmp = $chunks[1];
          if($i > $ndim && $i <= ($ndim+$nw)) { $vars[] = $var; }
        }
        $i++;
        if($i > $ndim && $i <= ($ndim+$nw)) { $vars[] = $tmp; }
        return $vars;
      }
    }
  }
  if($plottype == "2Didl--old") {
    Exec("strings $filename | head -3", $return);
    foreach ($return as $tmp) {
      if(eregi("rho", $tmp)) {
        $tmp = trim($tmp);
        $tmp = trim($tmp, 'O');
        $tmp = trim($tmp);
        while(eregi(" ", $tmp)) {
          $chunks = split(" ", $tmp, 2);
          $var = $chunks[0];
          $tmp = $chunks[1];
          $vars[] = $var;
        }
        $vars[] = $tmp;
        return $vars;
      }
    }
  }
  Exec("$pltview $filename", $return);
  foreach ($return as $tmp) {
    if(eregi("Var Names", $tmp)) {
      $tmp = trim($tmp, 'Var Names   : ');
      $tmp = trim($tmp);
	   while(eregi(",", $tmp)) {
        $chunks = split(",", $tmp, 2);
	     $var = $chunks[0];
	     $tmp = $chunks[1];
        $vars[] = $var;
      }
      $vars[] = $tmp;
      return $vars;
    }
  }
}

function PrintListMenu($list, $default, $name, $format) { // format = "name" or "number"
   echo "<SELECT NAME=$name>\n";
   // See if default exists
   $i = 0;
   $found = "0";
   foreach ($list as $value) {
      $i++;
      if($i == 1) {
         if($format == "name") {
            $newdefault = $value;
         } else {
            $newdefault = $i;
         }
      } else {
         if($format == "name") {
	    if($value == $default) { $newdefault = $default; }
         } else {
	    if($i == $default) { $newdefault = $default; }
         }
      }
   }
   // Write out menu
   $i = 0;
   foreach ($list as $value) {
      $i++;
      $option = $i;
      if($format == "name") { $option = $value; }
      echo "   <OPTION VALUE=\"$option\"";
      if ($option == $newdefault) {
        echo " SELECTED";
      }
      $prefix = "V$i: ";
      if($format == "name") { $prefix = ""; }
      echo ">$prefix$value</OPTION>\n";
   }
   echo '</SELECT>';
}

function PrintPlotMenu($list, $name, $prefix) {
   echo "<SELECT onchange=\"window.location=this.value\">\n";
   echo "   <OPTION VALUE=\"\"> -SELECT PLOT- </OPTION>\n";
   // Write out menu
   foreach ($list as $value) {
      echo "   <OPTION VALUE=\"$prefix$value\">$value</OPTION>\n";
   }
   echo '</SELECT>';
}

function PrintIsoVarsMenu($vars, $default) {
   echo "<SELECT NAME=isovar>\n";
   echo "   <OPTION VALUE=\"\">- REQUIRED IF YES &nbsp;</OPTION>\n";
   $i = 0;
   foreach ($vars as $value) {
      $i++;
      echo "   <OPTION VALUE=\"$i\"";
      if ($i == $default) {
        echo " SELECTED";
      }
      echo ">$value</OPTION>\n";
   }
   echo '</SELECT>';
}

function GetVectorVar($vector, $vars) {
   $value = 0;
   $i = 0;
   foreach ($vars as $value) {
      $i++;
		// tecplot
      if(eregi("${vector}_x", $value)) {
         $found = $i;
         return $found;
      }
		// idl
      if(eregi("${vector}x", $value)) {
         $found = $i;
         return $found;
      }
   }
   return NULL;
}

function SetLastPlotfile($plotfilelist, &$plotfile) {
// For a given plotfile, return the next and previous ones
   $plotfile = NULL;
   $lastplotfile = count($plotfilelist) - 1;
   if ($lastplotfile < 0) { return; }
   $plotfile = $plotfilelist[$lastplotfile];
   return;
}

function PrintComponents() {
   global $scriptname, $runname, $runpath;
   $common = "$scriptname?runname=$runname&cmp=";
   echo "
<table cellpadding=\"5\" cellspacing=\"0\" border=\"0\">
<tr>
  <td width=\"100\" class=\"cmp\">COMPONENT</td>
   ";
   if (is_dir("$runpath/SC")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}SC\" title=\"Solar Corona\">SC</a></td>"; }
   if (is_dir("$runpath/IH")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}IH\" title=\"Inner Heliosphere\">IH</a></td>"; }
   if (is_dir("$runpath/SP")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}SP\" title=\"Solar Energetic Particles\">SP</a></td>"; }
   if (is_dir("$runpath/GM")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}GM\" title=\"Global Magnetosphere\">GM</a></td>"; }
   if (is_dir("$runpath/IE")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}IE\" title=\"Ionospheric Electrodynamics\">IE</a></td>"; }
   if (is_dir("$runpath/IM")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}IM\" title=\"Inner Magnetosphere\">IM</a></td>"; }
   if (is_dir("$runpath/UA")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}UA\" title=\"Upper Atmosphere\">UA</a></td>"; }
   if (is_dir("$runpath/RB")) {
     echo "  <td width=\"30\" class=\"cmp\"><a href=\"${common}RB\" title=\"Radiation Belts\">RB</a></td>"; }
   echo "
</tr>
</table>
<br>
   ";
}

function PrintPltForm ($runname, $plotfilelist, $plotfile=0, $plotstyle) {
   global $scriptname, $cmp, $plottype;
   $common = "$scriptname?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=$plotstyle&plotfile=";
   echo "<SELECT onchange=\"window.location='$common'+this.value\">\n";
   foreach ($plotfilelist as $value) {
      echo "   <OPTION value=\"$value\"";
      if ($value == $plotfile) {
         echo " SELECTED";
      }
      $name = decodeFilename($value);
      echo ">$name</OPTION>\n";
   }
   echo "</SELECT>";
   echo "<br><br>\n";
}

function PrintXYZPltForm ($runname, $plotfile=0, $plotstyle) {
   global $scriptname, $cmp, $plottype, $runpath;
   $common = "$scriptname?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=$plotstyle&plotfile=";
   $choices = array('X', 'Y', 'Z', 'C2', 'C3');
   foreach ($choices as $choice) {
      $plotfilelist2 = GetXYZPlotList("$runpath/$cmp", "$choice");
      $countfiles2 = count($plotfilelist2);
      if($countfiles2 > 0) {
         echo "<SELECT onchange=\"window.location='$common'+this.value\">\n";
         echo "   <OPTION value=\"\" SELECTED>$choice: --not selected--</OPTION>\n";
         foreach ($plotfilelist2 as $value) {
            echo "   <OPTION value=\"$value\"";
            if ($value == $plotfile) {
               echo " SELECTED";
            }
            $name = decodeFilename($value);
            echo ">$choice: $name</OPTION>\n";
         }
         echo "</SELECT>";
      }
	}
   echo "<br><br>\n";
}

function decodeFilename($value) {
// Create more easily understood filename
   if(eregi(".plt", $value)) {
     // tecplot file.
     $chunks = split(".plt", $value, 2);
     $tmp = $chunks[0];
     if(eregi("sat", $tmp)) {
       // satelite file
       $chunks = split("sat_", $tmp, 2);
       $tmp = $chunks[1];
       $name = "$tmp";
       $name = trim($name);
       return $name;
     }
     $iter = "";
     if(eregi("_n", $tmp)) {
       $chunks = split("_n", $tmp, 2);
       $tmp = $chunks[0];
       $iter = "N=$chunks[1]";
     }
     $time = "";
     if(eregi("_t", $tmp)) {
       $chunks = split("_t", $tmp, 2);
       $time = $chunks[1];
       $timeH = substr($time,0,4);
       $timeM = substr($time,4,2);
       $timeS = substr($time,6,2);
       $time = "T=$timeH:$timeM:$timeS";
     }
     if($time || $iter) {
       $name = "$time $iter";
     } else {
       $name = "$value";
     }
     $name = trim($name);
     return $name;
   }
   if(eregi(".idl", $value)) {
     // idl IE file.
     $chunks = split(".idl", $value, 2);
     $tmp = $chunks[0];
     $chunks = split("it", $value, 2);
     $tmp = $chunks[1];
     $time = "";
     if(eregi("_", $tmp)) {
       $chunks = split("_", $tmp, 3);
       $ymd = $chunks[0];
       $hms = $chunks[1];
       $timeYe = substr($ymd,0,2);
       if($timeYe < 50) { $timeYe = "20$timeYe"; } else { $timeYe = "19$timeYe"; }
       $timeMo = substr($ymd,2,2);
       $timeDa = substr($ymd,4,2);
       $timeHo = substr($hms,0,2);
       $timeMi = substr($hms,2,2);
       $timeSe = substr($hms,4,2);
       $time = "T=$timeYe:$timeMo:$timeDa - $timeHo:$timeMi:$timeSe";
     }
     if($time) {
       $name = "$time";
     } else {
       $name = "$value";
     }
     $name = trim($name);
     return $name;
   }
   if(eregi(".out", $value)) {
     // idl file.
     $chunks = split(".out", $value, 2);
     $tmp = $chunks[0];
     $iter = "";
     if(eregi("_n", $tmp)) {
       $chunks = split("_n", $tmp, 2);
       $tmp = $chunks[0];
       $iter = "N=$chunks[1]";
     }
     $time = "";
     if(eregi("_t", $tmp)) {
       $chunks = split("_t", $tmp, 2);
       $time = $chunks[1];
       $timeH = substr($time,0,4);
       $timeM = substr($time,4,2);
       $timeS = substr($time,6,2);
       $time = "T=$timeH:$timeM:$timeS";
     }
     if($time || $iter) {
       $name = "$time $iter";
     } else {
       $name = "$value";
     }
     $name = trim($name);
     return $name;
   }
}

?>
