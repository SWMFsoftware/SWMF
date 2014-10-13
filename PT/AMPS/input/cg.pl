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


my $InputFileName="cg.input.Assembled.Block";   #$ARGV[0];
my $SpeciesFileName="cg.input.Assembled.Species";
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
my $spec;


#get the list of the species                                                   
my $TotalSpeciesNumber;
my @SpeciesList;

open (SPECIES,"<$SpeciesFileName") || die "Cannot open file $SpeciesFileName\\
n";
$TotalSpeciesNumber=<SPECIES>;
@SpeciesList=<SPECIES>;
close (SPECIES);

chomp($TotalSpeciesNumber);
foreach (@SpeciesList) {
    chomp($_);
}


my @BjornSourceRate=(0)x$TotalSpeciesNumber;
my @UniformSourceRate=(0)x$TotalSpeciesNumber;
my @JetSourceRate=(0)x$TotalSpeciesNumber;

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
  elsif ($InputLine eq "GRAVITY3D") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__3DGRAVITY__MODE_","_PIC_MODEL__3DGRAVITY__MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__3DGRAVITY__MODE_","_PIC_MODEL__3DGRAVITY__MODE__OFF_","pic/picGlobal.dfn");
      }
  }
  elsif ($InputLine eq "RADIATIVECOOLINGMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "CROVISIER") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIATIVECOOLING__MODE_","_PIC_MODEL__RADIATIVECOOLING__MODE__CROVISIER_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIATIVECOOLING__MODE_","_PIC_MODEL__RADIATIVECOOLING__MODE__OFF_","pic/picGlobal.dfn");
      }
  }
  elsif ($InputLine eq "RADIALVELOCITYMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIAL_VELOCITY_MODE_","_PIC_MODEL__RADIAL_VELOCITY_MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIAL_VELOCITY_MODE_","_PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_","pic/picGlobal.dfn");
      }
  }
  elsif ($InputLine eq "HELIOCENTRICDISTANCE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double HeliocentricDistance","$InputLine","main/Comet.cpp");
      ampsConfigLib::ChangeValueOfVariable("double HeliocentricDistance","$InputLine","main/main.cpp");
  }
  elsif ($InputLine eq "SUBSOLARPOINTAZIMUTH") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointAzimuth","$InputLine","main/Comet.cpp");
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointAzimuth","$InputLine","main/main.cpp");
  }
  elsif ($InputLine eq "SUBSOLARPOINTZENITH") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointZenith","$InputLine","main/Comet.cpp");
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointZenith","$InputLine","main/main.cpp");
  }
  elsif ($InputLine eq "NDIST") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("static int ndist","$InputLine","main/Comet.h");
  }
  elsif ($InputLine eq "BJORNPRODUCTIONRATEUSERDEFINED") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_","_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_","_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_OFF_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "BJORNPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $spec=$InputLine;
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $BjornSourceRate[getSpeciesNumber($spec)]=$InputLine;
      ampsConfigLib::ChangeValueOfArray("static double Bjorn_SourceRate\\[\\]",\@BjornSourceRate,"main/Comet.h");
  }
  elsif ($InputLine eq "UNIFORMPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $spec=$InputLine;
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $UniformSourceRate[getSpeciesNumber($spec)]=$InputLine;
      ampsConfigLib::ChangeValueOfArray("static double Uniform_SourceRate\\[\\]",\@UniformSourceRate,"main/Comet.h");
  }
  elsif ($InputLine eq "JETPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $spec=$InputLine;
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $JetSourceRate[getSpeciesNumber($spec)]=$InputLine;
      ampsConfigLib::ChangeValueOfArray("static double Jet_SourceRate\\[\\]",\@JetSourceRate,"main/Comet.h");
  }
  elsif ($InputLine eq "DUSTRMIN") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustSizeMin","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "DUSTRMAX") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustSizeMax","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "NDUSTRADIUSGROUPS") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("int DustSampleIntervals","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "DUSTTOTALMASSPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustTotalMassProductionRate","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "POWERLAWINDEX") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustSizeDistribution","$InputLine","main/Comet.cpp");
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
    my $species=$_[0];

    for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
	if ($species eq $SpeciesList[$i]) {
	    $res=$i;
	    last;
	}
    }

    if ($res eq -1) {
	die "Cannot fild species $_\n";
    }

    return $res;
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
