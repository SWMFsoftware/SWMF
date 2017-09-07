#!/usr/bin/perl



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
  

  if ($InputLine eq "TEST") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
     
    if ($InputLine eq "COLUMN_INTEGRATION_TEST") {
      ampsConfigLib::RedefineMacro("_TEST_","_COLUMN_INTEGRATION_TEST_","main/config.dfn");
    } 
    elsif ($InputLine eq "OFF") {
      ampsConfigLib::RedefineMacro("_TEST_","_OFF_","main/config.dfn");
    }
    else {
      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }

  #define gas mode or dust mode
  elsif ($InputLine eq "MODE") { 
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
      if ($InputLine eq "GAS") {
	  ampsConfigLib::RedefineMacro("_MODE_","_GAS_MODE_","main/config.dfn");
      } 
      elsif ($InputLine eq "DUST") {
	  ampsConfigLib::RedefineMacro("_MODE_","_DUST_MODE_","main/config.dfn");
      }
      else {
      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
      }
  }

    #load trajectory files or not
  elsif ($InputLine eq "LOADTRAJECTORYFILES") { 
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_LOAD_TRAJECTORY_FILES_","_ON_","main/main.cpp");
      } 
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_LOAD_TRAJECTORY_FILES_","_OFF_","main/main.cpp");
      }
      else {
      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
      }
  }

    #calculate surface exposure time or not
  elsif ($InputLine eq "SURFACEEXPOSURETIME") { 
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_SURFACE_EXPOSURE_TIME_CALCULATION_","_ON_","main/main.cpp");
      } 
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_SURFACE_EXPOSURE_TIME_CALCULATION_","_OFF_","main/main.cpp");
      }
      else {
      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
      }
  }



  elsif ($InputLine eq "OUT") {  #the output file number
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/[=();]/ /g;
      
      ($s0,$s1,$s2)=split(' ',$InputLine,3);
      ampsConfigLib::RedefineMacro("_OUTPUT_FILE_NUMBER_",$s1,"main/config.dfn");
  }
  
  elsif ($InputLine eq "NGAS") {  #the gas species number
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/[=();]/ /g;
      
      ($s0,$s1,$s2)=split(' ',$InputLine,3);
      ampsConfigLib::RedefineMacro("_GAS_SPEC_NUMBER_",$s1,"main/config.dfn");
  }
  
  elsif ($InputLine eq "NDUST") {
#      print "$InputLine, $InputComment \n";
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
#      print "$InputLine \n";
      my $ndust;
      $ndust=$InputLine;
      chomp($line);
      $line=~s/[();]/ /g;
      $line=~s/(=)/ /;
      $line=~s/(=)/ /;

      my $s0;
      ($s0,$line)=split(' ',$line,2);
      ($s0,$line)=split(' ',$line,2);

      while (defined $line) {
          ($s0,$line)=split(' ',$line,2);
          $s0=uc($s0);
          if ($s0 eq "NDUSTGROUP") {
	      ($s0,$line)=split(' ',$line,2);	  
	      ampsConfigLib::RedefineMacro("_DUST_SPEC_NUMBER_",$ndust,"main/config.dfn");
	      ampsConfigLib::RedefineMacro("_DUST_GROUP_NUMBER_",$s0,"main/config.dfn");
	      ampsConfigLib::RedefineMacro("_DUST_CASE_","_DUST_CASE__".$ndust."SPEC_".$s0."GROUP_","main/config.dfn");
 
          }
          else {
	      warn ("Cannot recognize the option (line=$InputLine, nline=$InputFileLineNumber)");
	      die "$InputLine: Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
          }
      }
      
  }

#   elsif ($InputLine eq "SPICEKERNELSPATH") {
#       ($InputLine,$InputComment)=split('!',$line,2);
#       chomp($InputLine);
#       $InputLine=~s/=/ /g;
      
#       ($InputLine,$InputComment)=split(' ',$InputLine,2);
#       chomp($InputComment);
      
#       $InputComment="('".$InputComment."')";
#       ampsConfigLib::ChangeValueOfVariableNoSemiColon("PATH_VALUES",$InputComment,"main/kernels.tm");  
#   }
 
   
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
  
}  
 

my $file = '.amps.conf'; 
open my $info,$file || die "Cannot find file \".amps.conf\"\n";

while ($line=<$info>) {

  $LineOriginal=$line;
  chomp($line);
  
  
  ($InputLine,$InputComment)=split('!',$line,2);
  $InputLine=uc($InputLine);
  chomp($InputLine);
 
  $InputLine=~s/[=()]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;

  if ($InputLine eq "SPICEKERNELS") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/=/ /g;
      
      ($InputLine,$InputComment)=split(' ',$InputLine,2);
      chomp($InputComment);
      
      $InputComment="('".$InputComment."')";
      ampsConfigLib::ChangeValueOfVariableNoSemiColon("PATH_VALUES",$InputComment,"main/kernels.tm");  
  }


}


