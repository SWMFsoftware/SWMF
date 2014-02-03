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


print "Process the exospehre model input:\n";


my $InputFileName="mercury.input.Assembled.Block";   #$ARGV[0];
my $SpeciesFileName="mercury.input.Assembled.Species";
my $WorkingSourceDirectory="src";  #$ARGV[1];

$ampsConfigLib::WorkingSourceDirectory=$WorkingSourceDirectory;

my $line;
my $LineOriginal;
my $FileName;

my $InputFileLineNumber=0;
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
    my $HeliocentricDistance="_AU_";
    my $AzimuthCenter;
    my $ZenithCenter;
    my $Angle;
    my @JetSourceRate=(0)x$TotalSpeciesNumber;
    my @BjornSourceRate=(0)x$TotalSpeciesNumber;
    my $spec;
    my $ndist;
    my $ct=0;

    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    $LineOriginal=$line;
    $InputFileLineNumber++;
    chomp($line);
    
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    
    
    $InputLine=~s/[=()]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);
    $InputLine=~s/ //g;
    
    while (defined $InputComment) {
	if ($ct>0) {
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	}
	$ct=$ct+1;
	if ($InputLine eq "HELIOCENTRICDISTANCE") {
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $HeliocentricDistance=$InputLine    ;
	    printf("HELIOCENTRICDISTANCE = ");
	    printf("$InputLine \n");
	    ampsConfigLib::ChangeValueOfVariable("static const double ImpactVaporization_HeliocentricDistance",$HeliocentricDistance,"models/exosphere/Exosphere.h");
	}
	elsif ($InputLine eq "AZIMUTHCENTER") { 
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $AzimuthCenter=$InputLine;    
	    printf("AZIMUTHCENTER = ");
	    printf("$InputLine \n");
	    ampsConfigLib::ChangeValueOfVariable("static double azimuthCenter",$AzimuthCenter,"main/Comet.cpp");
	}
	elsif ($InputLine eq "ZENITHCENTER") { 
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $ZenithCenter=$InputLine;  
	    printf("ZENITHCENTER = ");
	    printf("$InputLine \n");
	    ampsConfigLib::ChangeValueOfVariable("static double zenithCenter",$ZenithCenter,"main/Comet.cpp");
	}
	elsif ($InputLine eq "ANGLE") { 
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $Angle=$InputLine;
	    printf("ANGLE = ");
	    printf("$InputLine \n");
	    ampsConfigLib::ChangeValueOfVariable("static double angle",$Angle,"main/Comet.cpp");
	}
	elsif ($InputLine eq "JETPRODUCTIONRATE") {
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $spec=$InputLine;
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $JetSourceRate[getSpeciesNumber($spec)]=$InputLine;
	    printf("JETPRODUCTIONRATE = ");
	    printf("$InputLine \n");
	    ampsConfigLib::ChangeValueOfArray("static double Jet_SourceRate\\[\\]",\@JetSourceRate,"main/Comet.h");
	}
	elsif ($InputLine eq "NDIST") { 
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $ndist=$InputLine;    
	    printf("NDIST = ");
	    printf("$InputLine \n");
	    ampsConfigLib::ChangeValueOfVariable("static int ndist",$ndist,"main/Comet.h");
	}
	elsif ($InputLine eq "BJORNPRODUCTIONRATE") {
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $spec=$InputLine;
	    ($InputLine,$InputComment)=split(' ',$InputComment,2);
	    $InputLine=~s/ //g;
	    $BjornSourceRate[getSpeciesNumber($spec)]=$InputLine;
	    printf("BJORNPRODUCTIONRATE = ");
	    printf("$InputLine \n");
	    ampsConfigLib::ChangeValueOfArray("static double Bjorn_SourceRate\\[\\]",\@BjornSourceRate,"main/Comet.h");
	}
	else {
	    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
	}
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
    die "Cannot fild species $_\n";
  }
  
  return $res;
}


=comment
#=============================== Change a value of a variable in the code  =============================
sub ChangeValueOfVariable {
  my $Variable=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$WorkingSourceDirectory/$File") || die "Cannot open file $WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$WorkingSourceDirectory/$File");

  foreach (@FileContent) {  
    if ($_=~m/($Variable)/) {
      my $t=$Variable;
      
      $t=~s/\\//g;
      $_=$t."=".$Value.";\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

#=============================== Change a value of a array in the code  =============================
sub ChangeValueOfArray {
  my $Variable=$_[0];
  my @Value=@{$_[1]};
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$WorkingSourceDirectory/$File") || die "Cannot open file $WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    if ($_=~m/($Variable)/) {
      my $t=$Variable;
      
      $t=~s/\\//g;
      $_=$t."={";
      
      for (my $i=0;$i<1+$#Value;$i++) {
        $_=$_."$Value[$i]";
        
        if ($i ne $#Value) {
          $_=$_.",";
        }
      }
      
      $_=$_."};\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT); 
}

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
      print "\#define $Macro $Value\n";
      $_="\#define $Macro $Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}
=cut
