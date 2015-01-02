#!/usr/bin/perl
#$Id$
#reads the "ccmc" section of the input file 

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


my $InputFileName=$ARGV[0];
my $WorkingSourceDirectory=$ARGV[1];
my $InputDirectory=$ARGV[2];  #the location of the data files that will be used in the model run

$ampsConfigLib::WorkingSourceDirectory=$WorkingSourceDirectory;

my $line;
my $LineOriginal;
my $InputFileLineNumber=0;
my $FileName;
my $InputLine;
my $InputComment;


open (InputFile,"<","$InputFileName") || die "Cannot find file \"$InputFileName\"\n";

while ($line=<InputFile>) {
  ($InputFileLineNumber,$FileName)=split(' ',$line);
  $line=<InputFile>;
  
  $LineOriginal=$line;
  chomp($line);
  
  
  ($InputLine,$InputComment)=split('!',$line,2);
  $InputLine=uc($InputLine);
  chomp($InputLine);
 
  $InputLine=~s/[=(),]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  
  if ($InputLine eq "XMIN") {
    my $i;
    my @xmin=(0)x3;
    
    for ($i=0;$i<3;$i++) {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $xmin[$i]=$InputLine;
    }
    
    ampsConfigLib::ChangeValueOfArray("static const double xmin\\[\\]",\@xmin,"main/ccmc.h");
  }
  elsif ($InputLine eq "XMAX") {
    my $i;
    my @xmax=(0)x3;
    
    for ($i=0;$i<3;$i++) {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $xmax[$i]=$InputLine;
    }
    
    ampsConfigLib::ChangeValueOfArray("static const double xmax\\[\\]",\@xmax,"main/ccmc.h");    
  }
  elsif ($InputLine eq "SPHERE") {
    while (defined $InputComment) {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
          
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CCMC_INTERNAL_BOUNDARY_MODE_","_PIC_MODE_ON_","main/ccmc.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CCMC_INTERNAL_BOUNDARY_MODE_","_PIC_MODE_OFF_","main/ccmc.dfn")
      }
      elsif ($InputLine eq "RADIUS") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        ampsConfigLib::AddLine2File("#undef _CCMC_INTERNAL_BOUNDARY_RADIUS_\n#define _CCMC_INTERNAL_BOUNDARY_RADIUS_ $InputLine\n\n","main/ccmc.dfn")
      }
      else {
        die "Option is unknown #1 ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
  }
  elsif ($InputLine eq "#ENDBLOCK") {
      last;
  }
  else {
    die "Option is unknown #2 ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
  }
}


