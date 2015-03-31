#!/usr/bin/perl
#$Id$
#reads the "ChargeExchange" section of the input file 

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
my ($Surface,$Model);

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
  
  if ($InputLine eq "MODEL") {
    ($Model,$InputComment)=split(' ',$InputComment,2);    
  }
  elsif ($InputLine eq "#ENDBLOCK") {
    #update the model settings
    my $s;
    
    if ((!defined $Model)) {
      print "Some of the model parameters are not defined!!!!!! \n";
      die;
    }
    
    ampsConfigLib::AddLine2File("#undef _CHARGE_EXCHANGE__MODEL_\n#define _CHARGE_EXCHANGE__MODEL_ _CHARGE_EXCHANGE__".$Model."_","models/charge_exchange/ChargeExchange.dfn");
    
    #quit the loop
    last;
  }
  else {
    die "Option is unknown #2 ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
  }
}


