#!/usr/bin/perl
#$Id$
#average multiple output files produced by AMPS

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);

#$nNode,$nCells -> the number of the nodes and cells, $nTitleLines -> the number of the title lines (TITLE, ZODE, VARABLE) in the output file

my ($str,$s0,$s1,$s2,$line,$nNode,$nCells,$nTitleLines,$i,$nDataFiles,$nVariables);
my (@data,@sum);

#arguments:
#$ARGV[0] -> the name of the intput file
#$ARGV[1] -> the working directory
print "Average AMPS output files:\nArguments: @ARGV\n";


#the output file
open (fAveraged,">averaged.dat");

#read the variable list and the number of the number of the cells and nodes

open (fInput,"<$ARGV[0]");
$nTitleLines=0;

for ($i=0;$i<3;$i++) {
  $line=<fInput>;
  chomp($line);
  $str=uc($line);  
  $str=~s/[(),=]/ /g;

  ($s0,$s1)=split(' ',$str,2); 
 
  if ($s0 eq "VARIABLES") {
    $nTitleLines++;
    print fAveraged "$line \n";
  }
  elsif ($s0 eq "TITLE") {
    $nTitleLines++;
    print fAveraged "$line\n";
  }
  elsif ($s0 eq "ZONE") {
    $nTitleLines++;
    print fAveraged "$line\n";

    while (defined $s1) {
      ($s0,$s2,$s1)=split(' ',$s1,3);
 
      if (defined $s0) {
        if ($s0 eq "E") {
          $nCells=$s2;
        }
        elsif ($s0 eq "N") {
          $nNode=$s2;
        }
      }
    }
  }
}    

#get the number of the variables
$line=<fInput>;
@data=split(' ',$line);
$nVariables=scalar(@data);

print "Number of cells: $nCells\n";
print "Number of nodes: $nNode\n";
print "Number of variables: $nVariables\n";

close (fInput);

#open the input data files
my @filehandles;

print "Files to avarage: @ARGV \n";

$nDataFiles=scalar(@ARGV);
print "Output file number: $nDataFiles\n";

for ($i=0;$i<$nDataFiles;$i++) {
  print "$ARGV[$i]\n";

  open(FILE, "<$ARGV[$i]") || die;
  #push the typeglobe to the end of the array
  push(@filehandles, *FILE);    

  for ($i=0;$i<$nTitleLines;$i++) {
    $line=<FILE>;
  }
}



#average of output the model data
my ($nfile,$nline,$file);
my @tmp=((0)x$nVariables);

for ($nline=0;$nline<$nNode;$nline++) {   
  for ($i=0;$i<$nVariables;$i++) {
    $tmp[$i]=0.0;
  }

  foreach $file (@filehandles) {
    $line=<$file>;
    @data=split(' ',$line);

    for ($i=0;$i<$nVariables;$i++) {
      $tmp[$i]+=$data[$i];
    }
  }

  for ($i=0;$i<$nVariables;$i++) {
    $tmp[$i]/=$nVariables;
  }

  print fAveraged "@tmp\n";
} 

print "#1\n";

#output the connectivity list
$file=$filehandles[0];

for ($nline=0;$nline<$nCells;$nline++) {
  $line=<$file>;
  print fAveraged "$line\n";
}



#close open file
foreach $file (@filehandles) {
  close ($file);
}

close (fAveraged);
   






