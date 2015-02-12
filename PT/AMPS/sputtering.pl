#!/usr/bin/perl
#reads the "sputtering" section of the input file 

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


#print "Process the sputtering model input:\n";


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
  
  if ($InputLine eq "SPUTTERINGMODE") {
      while (defined $InputComment) {
	  ($InputLine,$InputComment)=split(' ',$InputComment,2);
          
	  if ($InputLine eq "ON") {
	      ampsConfigLib::RedefineMacro("_SPUTTERING__MODE_","_SPUTTERING__ON_","models/sputtering/Sputtering.dfn");
	  }
	  elsif ($InputLine eq "OFF") {
	      ampsConfigLib::RedefineMacro("_SPUTTERING__MODE_","_SPUTTERING__OFF_","models/sputtering/Sputtering.dfn");
	  }
	  else {
	      die "Option is unknown #1 ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
	  }
      }
  }
  elsif ($InputLine eq "SPUTTERINGSURFACE") {    
      while (defined $InputComment) {
	  ($InputLine,$InputComment)=split(' ',$InputComment,2);
          
	  if ($InputLine eq "ICE") {
	      ampsConfigLib::RedefineMacro("_SPUTTERING__SURFACE_","_SPUTTERING__ICE_","models/sputtering/Sputtering.dfn");
	      ($InputLine,$InputComment)=split(' ',$InputComment,2);
	      if ($InputLine eq "MODEL") {
		  ($InputLine,$InputComment)=split(' ',$InputComment,2);
		  if ($InputLine eq "RUBIN") {
		      ampsConfigLib::RedefineMacro("_SPUTTERING__ICE_MODEL_","_SPUTTERING__ICE__RUBIN_","models/sputtering/Sputtering.dfn");
		  }
		  elsif ($InputLine eq "TEOLIS") {
		      ampsConfigLib::RedefineMacro("_SPUTTERING__ICE_MODEL_","_SPUTTERING__ICE__TEOLIS_","models/sputtering/Sputtering.dfn");
		  }
		  else {
		      die "Option is unknown #1 ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
		  }
	      }
	      else {
		  die "Option is unknown #1 ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
	      }
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
