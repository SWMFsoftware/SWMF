<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php    // Create a list of thumbnails for plot style

   if(! is_dir("$runpath/images")) {
     Exec("cd $runpath; mkdir images");
   }
   if ($cmp) {
     if(! $plottype) { // Select plottype before style
       Exec("rsync -a --exclude \"CVS\" BASE/$cmp* $runpath/images/");
       $types = array();
       $tmpdir = opendir( "$runpath/images" );
       while( $type = readdir( $tmpdir ) ) {
         if (ereg("${cmp}_", $type)) { $types[] = $type; }
       }
       if($types){ sort($types);
         echo "<center><H2>Select $cmp Filetype</H2></center>";
         echo "<ul>\n";
         foreach ($types as $type) {
           $explode = explode("_", $type, 2); $name = $explode[1];
           $plottype = $name;
           $plotfilelist = GetPlotList("$runpath/$cmp");
           $countfiles = count($plotfilelist);
           if($countfiles < 1) {
             echo "<li><b>$name</b> (NO files)</li>\n";
           } else {
             echo "<li><a href=\"$scriptname?runname=$runname&cmp=$cmp&plottype=$name&plotstyle=001&plotfile=$plotfile\"><b>$name</b></a> ($countfiles files)</li>\n";
           }
           $plottype = NULL;
         }
         echo "</ul>\n";
       }
     }
     if($plottype) {
       echo "<center><H3>Plottype: $plottype</H3></center>";
       $plotstylelist = GetStyleList("$runpath/images/${cmp}_$plottype");
       echo "<center><H2>$cmp Plot Styles</H2></center>";
       echo "<P ALIGN=\"RIGHT\">";
       foreach ($plotstylelist as $style) {
         $pieces = explode(".", $style);
         $plotstyleclip = $pieces[0];
         $pieces = explode("-", $plotstyleclip);
         $number = $pieces[1];
         if ($number == $plotstyle) { 
           echo "<div class=\"selected\">";
         } else {
           echo "<div class=\"unselected\">";
         }
         //echo "$number\n";
         echo "<table><td width=35 align=left>$number</td>";
         echo "<td width=35 align=left><a href=\"movie.php?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=${number}&wait=1\" TARGET=\"_movie_$runname\">movie</a></td>";
         echo "  <td width=50 align=right><a href=\"delete.php?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=$number&plotfile=$plotfile\">delete</a></td></table>\n";
         $file = "../images/noplot.png";
         $dir = opendir("$runpath/images/${cmp}_$plottype");
         while( $newfile = readdir( $dir ) ) {
           if (eregi("-$number.png", $newfile)) {
             $file = "$runpath/images/${cmp}_$plottype/$newfile";
             break;
           }
         }
         echo "<A HREF=\"$scriptname?runname=$runname&cmp=$cmp&plottype=$plottype&plotstyle=$number&plotfile=$plotfile\">"; 
         echo "<IMG width=130 SRC=\"$file\" BORDER=0>";
         echo "</A><BR>";
         echo "</div>\n";
       }
       echo "</P>\n";
     }
   }

?>
