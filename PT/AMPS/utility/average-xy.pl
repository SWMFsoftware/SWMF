#!/usr/bin/perl
#$Id$
#average multiple output files produced by AMPS

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);
use IO::File;
use Scalar::Util qw(looks_like_number);

print "Average AMPS output files:\n";

#the output file
open (fAveraged,">amps.averaged.dat");

#get the number of variables in a line of the averaged file 
my ($i,$nVariables,$nDataFiles,$line,$s0,$s1);
my $flag=0;

open (fInput,"<$ARGV[0]");  

while ($line=<fInput>) {
  ($s0,$line)=split(' ',$line,2);
  
  if (looks_like_number($s0)) {
    #a data line is found, determine how name elements does it has
    $nVariables=1;
    chomp($line);

    while (defined $line) {
      ($s0,$line)=split(' ',$line,2);
      $nVariables+=1;  
    }

    $flag=1;
    last;
  }
}

close (fInput);

#check wether a data line has been found
if ($flag == 0) {
  die "No data lines was found in file $ARGV[0]\n";
}

#average the data 
$nDataFiles=scalar(@ARGV);

print "The number of the variables: $nVariables\n";
print "The number of the data files: $nDataFiles\n";

#create an array of the file handelers
my (@filehandles,$file);

for ($i=0;$i<$nDataFiles;$i++) {
  print "$ARGV[$i]\n";

  $file=IO::File->new("< $ARGV[$i]") || die "Cannot open file $ARGV[$i]\n"; 
  push(@filehandles,$file);
}

#avarage the data files 
my @tmp=((0)x$nVariables); 
my ($iVar,$iFile,$LineOriginal);
my $FistDataFile;
my $OtherDataFile;

$FistDataFile=$filehandles[0];

while ($line=<$FistDataFile>) {
  #check whether the line contains numbers
  $s1=$line;
  $LineOriginal=$line;
  ($s0,$s1)=split(' ',$s1,2);  

  if (looks_like_number($s0)) {
    #the line contains numbers -> average them 
    for ($iVar=0;$iVar<$nVariables;$iVar++) {
      ($s0,$line)=split(' ',$line,2);
      $tmp[$iVar]=$s0;
    }

    for ($iFile=1;$iFile<$nDataFiles;$iFile++) {
      $OtherDataFile=$filehandles[$iFile];
      $line=<$OtherDataFile>;

      for ($iVar=0;$iVar<$nVariables;$iVar++) {
        ($s0,$line)=split(' ',$line,2);

        $tmp[$iVar]+=$s0;
       }  
    }

    #output the averaged data into a file
    for ($iVar=0;$iVar<$nVariables;$iVar++) {
      my $t;

      $t=$tmp[$iVar]/$nDataFiles;
      print fAveraged  "$t ";
    } 

    print fAveraged  "\n";
  }
  else {
    #the line has no numbers -> just read the another line from the rest of the files 
    for ($iFile=1;$iFile<$nDataFiles;$iFile++) {
      $OtherDataFile=$filehandles[$iFile];
      $line=<$OtherDataFile>;
    }

    print fAveraged  "$LineOriginal\n";
  }
}

#close the output file 
close (fAveraged);
