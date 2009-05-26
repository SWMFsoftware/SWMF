#!/usr/bin/perl -w
# Script: replace_amrinit.pl
# Replaces "#AMRINIT" section in PARAM.in file by #USERINPUTBEGIN and
# saves initial PARAM.in file as PARAM.in_old_version
# One input parameter: Grid resolution
# if not specified, will use default value:0.03125
#
if ( ! -e "PARAM.in" ){
     print "Please run this script in the directory containing  
PARAM.in file\n";
     print  "exiting script\n";    exit;
}

# Determine which line contains AMRINIT identifier:
$line_with_amrinit=`grep -n "#AMRINIT" PARAM.in |head -1 |cut -f1 -d:`;

if ( $line_with_amrinit lt 0){
     print "PARAM.in file does not contain #AMRINIT block\n";
     print  "exiting script\n";
exit;
}

$line_with_amrinit= $line_with_amrinit - 1;
$numArgs = $#ARGV;
if ("$numArgs" lt "0" ){
     $gridresolution=0.03125;  # set default
     print "Grid resolution set to default value: $gridresolution\n";
} else {
     $gridresolution=$ARGV[0];
     print "Grid resolution requested: $gridresolution\n";
}


$line_with_grid_name= $line_with_amrinit+1;
$line_with_grid_level= $line_with_amrinit+2;

$first_line_after_amrinit_block=$line_with_amrinit + 3;
$last_line_in_file=`wc -l PARAM.in | awk '{print \$1}'`;
$final_lines = $last_line_in_file - $first_line_after_amrinit_block;

open (PARAMFILE,"<PARAM.in");
@logmes= <PARAMFILE>;
close(PARAMFILE);
@PARAMFILE=@logmes;

chomp($grid_name = `echo "$PARAMFILE[$line_with_grid_name]" | awk  
'{print \$1}'`);
chomp($grid_name = $grid_name);
chomp($grid_level = $PARAMFILE[$line_with_grid_level]);

$replacement_lines="#USERINPUTBEGIN

#CCMCGRID
$grid_name    InitialRefinementType (character)

#USERINPUTEND

#GRIDLEVEL
$grid_level
init

#GRIDRESOLUTION
$gridresolution
user
";

system("cp -p PARAM.in PARAM.in_old_version");

system("head -$line_with_amrinit PARAM.in_old_version > PARAM.in");

open (PARAMFILE,">>PARAM.in");
print PARAMFILE $replacement_lines;
close PARAMFILE;

system("tail -$final_lines PARAM.in_old_version >> PARAM.in");

exit;
