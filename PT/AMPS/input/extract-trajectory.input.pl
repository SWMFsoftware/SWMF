#!/usr/bin/perl
#$Id$
#reads the corresponding section of the input file and modifies source code

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);

use ampsConfigLib;


#my $command= $];
#print "Perl version : ".$command;



#arguments: 
#$ARGV[0] -> the name of the intput file 
#$ARGV[1] -> the working directory


#print "Process the exospehre model input:\n";


my $InputFileName=$ARGV[0]; #"europa.input.Assembled.Block"; #$ARGV[0];  # moon.input.Assembled.Block";
my $SpeciesFileName=$InputFileName; $SpeciesFileName =~ s/\.Block$/.Species/;
my $WorkingSourceDirectory=$ARGV[1];   #"srcTemp"; #$ARGV[1];   # srcTemp

$ampsConfigLib::WorkingSourceDirectory=$WorkingSourceDirectory;

my $line;
my $LineOriginal;

my $InputFileLineNumber=0;
my $FileName;
my $InputLine;
my $InputComment;
my $s0;
my $s1;
my $s2;


#get the list of the speces 
my $TotalSpeciesNumber;
my @SpeciesList;

open (SPECIES,"<$SpeciesFileName") || die "Cannot open file $SpeciesFileName\n";
$TotalSpeciesNumber=<SPECIES>;
@SpeciesList=<SPECIES>;
close (SPECIES);
      
chomp($TotalSpeciesNumber);
foreach (@SpeciesList) {
  chomp($_);
}


#add the model header to the pic.h
#open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
#print PIC_H "#include \"Europa.h\"\n";
#close (PIC_H);

open (InputFile,"<","$InputFileName") || die "Cannot find file \"$InputFileName\"\n";

while ($line=<InputFile>) {
  ($InputFileLineNumber,$FileName)=split(' ',$line);
  $line=<InputFile>;
  
  $LineOriginal=$line;
  chomp($line);
  
  
  ($InputLine,$InputComment)=split('!',$line,2);
  $InputLine=uc($InputLine);
  chomp($InputLine);
 
  $InputLine=~s/[=()]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
 
  if ($InputLine eq "FILE") {
      my $fname;
      
      $line=~s/[();,]/ /g;
      $line=~s/[=]/ /;
      chomp($line);
      
      ($InputLine,$fname,$InputComment)=split(' ',$line,2);
      chomp($fname);
      
      
      print "$fname\n";
      print "$line\n";

      ampsConfigLib::ChangeValueOfVariable("char fname\\[\\]","\"".$fname."\"","main/main.cpp");  
  }

  elsif ($InputLine eq "TRAJECTORYPOINTBUFFERLENGTH") { 
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=();]/ /g;
    
    ($s0,$s1,$s2)=split(' ',$InputLine,3);
    ampsConfigLib::ChangeValueOfVariable("int TrajectoryPointBufferLength",$s1,"main/main.cpp");   
  }
  
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
  
}


#=============================== Determine the species number  =============================
sub getSpeciesNumber {
  my $res=-1;
  my $spec=$_[0];
  
  for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
    if ($spec eq $SpeciesList[$i]) {
      $res=$i;
      last;
    }
  } 
  
  if ($res eq -1) {
    die "Cannot find species $spec\n";
  }
  
  return $res;
}

