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


my $InputFileName="hartley2.input.Assembled.Block";   #$ARGV[0];
my $SpeciesFileName="hartley2.input.Assembled.Species";
my $WorkingSourceDirectory="srcTemp";  #$ARGV[1];

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


#add the model header to the pic.h
open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
print PIC_H "#include \"Comet.h\"\n";
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
  
  if ($InputLine eq "SODIUMSTICKINGPROBABILITY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
     
    if ($InputLine eq "CONST") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_","_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_","main/Comet.cpp");
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY__CONSTANT_VALUE_","$InputLine","main/Comet.cpp");
    } 
    elsif ($InputLine eq "YAKSHINSKIY2005SS") {
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_","_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__YAKSHINSKY2005SS_","main/Comet.cpp");
    }
    else {
      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  elsif ($InputLine eq "SODIUMREEMISSIONFRACTION") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    
    ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY__REEMISSION_FRACTION_","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }
   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
}


=comment
#=============================== Change a definition of a macro in a source code  =============================
sub RedefineMacro {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$WorkingSourceDirectory/$File") || die "Cannot open file $WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$WorkingSourceDirectory/$File") || die "Cannot open file $WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro /) {
#      print "\#define $Macro $Value\n";
      $_="\#define $Macro $Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}
=cut
