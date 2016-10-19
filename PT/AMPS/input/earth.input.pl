#!/usr/bin/perl
#$Id$
#reads the "Exosphere" section of the input file and modify the exosphere model of the PIC code

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
open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
close (PIC_H);

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
  
  if ($InputLine eq "MESHSIGNATURE") {
    $line=~s/[=():]/ /g;
    ($InputLine,$line)=split(' ',$line,2);
    ($InputLine,$line)=split(' ',$line,2);

    ampsConfigLib::ChangeValueOfVariable("char Earth::Mesh::sign\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputLine."\"","main/Earth.cpp");
  }
  elsif ($InputLine eq "SEP") {
    ($s0,$s1)=split(' ',$InputComment,2);

    if ($s0 eq "ON") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_SEP__MODE_","_PIC_MODE_ON_","main/Earth.dfn");
    } 
    elsif($s0 eq "OFF") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_SEP__MODE_","_PIC_MODE_OFF_","main/Earth.dfn");
    }
    else {
      die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
    }    
  }
  elsif ($InputLine eq "GCR") {
    ($s0,$s1)=split(' ',$InputComment,2);

    if ($s0 eq "ON") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_GCR__MODE_","_PIC_MODE_ON_","main/Earth.dfn");
    } 
    elsif($s0 eq "OFF") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_GCR__MODE_","_PIC_MODE_OFF_","main/Earth.dfn");
    }
    else {
      die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
    }     
  }
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
  
}
