#!/usr/bin/perl
#reads the "ChemicalSolver" section of the input file 

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

#parameters of the model settings
my ($Source,$SCoeff,$Product,$PCoeff,$Rate,$Misc,$LHS,$RHS);

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
  
  if ($InputLine eq "REACTIONS") {
      # Read reaction equations (all written in one line)
      # the format of an equation is the following
      # Left-Hand Side (LHS) -> Right-Hand Side (RHS) : Rate
      # LHS/RHS: SCoeff1 * Source1 + ... / PCoeff * Product1 + ...
      while(($LHS,$InputComment)=split('->',$InputComment,2)){
	  ($SCoeff,$LHS)=split('\*',$LHS,2); $SCoeff=~tr/ //d;
	  ($Source,$LHS)=split('\+',$LHS,2); $Source=~tr/ //d;
	  $Source = '_'.$Source.'_SPEC_';
	  while($LHS){
	      ($Misc,$LHS)=split('\*',$LHS,2);$Misc=~tr/ //d;
	      $SCoeff=$SCoeff.','.$Misc;
	      ($Misc,$LHS)=split('\+',$LHS,2);$Misc=~tr/ //d;
	      $Source=$Source.',_'.$Misc.'_SPEC_';
	  }
	  ($RHS,$InputComment)=split(':',$InputComment,2);
	  ($PCoeff,$RHS)=split('\*',$RHS,2);$PCoeff=~tr/ //d;
	  ($Product,$RHS)=split('\+',$RHS,2);$Product=~tr/ //d;
	  $Product='_'.$Product.'_SPEC_';
	  while($RHS){
	      ($Misc,$RHS)=split('\*',$RHS,2);$Misc=~tr/ //d;
	      $PCoeff=$PCoeff.','.$Misc;
	      ($Misc,$RHS)=split('\+',$RHS,2);$Misc=~tr/ //d;
	      $Product=$Product.',_'.$Misc.'_SPEC_';
	  }
	  ($Rate,$InputComment)=split(' ',$InputComment,2);$Rate=~tr/ //d;
	  # SCoeff, Source, PCoeff, Product, Rate NEED TO BE WRITTEN 
	  # TO SOURCE CODE FILE
      }
  }
  elsif ($InputLine eq "#ENDBLOCK") {
    #update the model settings
    my $s;    
   
    #quit the loop
    last;
  }
  else {
    die "Option is unknown #2 ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
  }
}


