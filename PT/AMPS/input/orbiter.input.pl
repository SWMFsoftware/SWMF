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
 
  $InputLine=~s/[=(),]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  
  if ($InputLine eq "MESHSIGNATURE") {
    $line=~s/[=():]/ /g;
    ($InputLine,$line)=split(' ',$line,2);
    ($InputLine,$line)=split(' ',$line,2);

    ampsConfigLib::ChangeValueOfVariable("char Orbiter::Mesh::sign\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputLine."\"","main/Orbiter.cpp");
  }
  elsif ($InputLine eq "SURFACEMODEL") {
    $line=~s/[=():]/ /g;
    ($InputLine,$line)=split(' ',$line,2);
    ($InputLine,$line)=split(' ',$line,2);

    ampsConfigLib::ChangeValueOfVariable("char Orbiter::SurfaceModel::MeshFileName\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputLine."\"","main/Orbiter.cpp");
  }   
  elsif ($InputLine eq "SURFACEMESHSCALINGFACTOR") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    ampsConfigLib::ChangeValueOfVariable("double Orbiter::SurfaceModel::ScalingFactor",$InputLine,"main/Orbiter.cpp");         
  }
  elsif ($InputLine eq "SURFACEMODELFORMAT") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "CEA") {
      ampsConfigLib::ChangeValueOfVariable("int Orbiter::SurfaceModel::MeshFileFormat","Orbiter::SurfaceModel::MeshFileFormat_CEA","main/Orbiter.cpp");      
    }
    elsif ($InputLine eq "NASTRAN") {
      ampsConfigLib::ChangeValueOfVariable("int Orbiter::SurfaceModel::MeshFileFormat","Orbiter::SurfaceModel::MeshFileFormat_NASTRAN","main/Orbiter.cpp");           
    }
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  
  elsif ($InputLine eq "UPSTREAMVELOCITY") {
    my @v=(0)x3;
    
    ($v[0],$v[1],$v[2],$InputLine)=split(' ',$InputComment,4);
    ampsConfigLib::ChangeValueOfArray("double Orbiter::UpstreamBC::Velocity\\[3\\]",\@v,"main/Orbiter.cpp");
  }
  elsif ($InputLine eq "UPSTREAMNUMBERDENSITY") {
    my @NumberDensityTable=((0)x$TotalSpeciesNumber);
    
    while (defined $InputComment) {
      my ($c0,$c1,$c2,$nspec);
      
      ($c0,$c1,$c2,$InputComment)=split(' ',$InputComment,4);
      
      if ($c0 eq "CONST") {
        $nspec=ampsConfigLib::GetElementNumber($c1,\@SpeciesList);
        
        if ($nspec==-1) {
          die "Cannot recognize species '$c1' in line $InputFileName ($line)\n";
        }
        
        $NumberDensityTable[$nspec]=$c2;
      }
    }

    ampsConfigLib::ChangeValueOfArray("double Orbiter::UpstreamBC::NumberDensity\\[PIC::nTotalSpecies\\]",\@NumberDensityTable,"main/Orbiter.cpp");
  }
  elsif ($InputLine eq "UPSTREAMTEMPERATURE") {
    my $t;
    
    ($t,$InputLine)=split(' ',$InputComment,2);
    ampsConfigLib::ChangeValueOfVariable("double Orbiter::UpstreamBC::Temeprature",$t,"main/Orbiter.cpp");
  } 
  elsif ($InputLine eq "DOMAINLIMITS") {
    my @xMinOffset=(0.5)x3;
    my @xMaxOffset=(0.5)x3;

    while (defined $InputComment) {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "XMINOFFSET") {
        ($xMinOffset[0],$xMinOffset[1],$xMinOffset[2],$InputComment)=split(' ',$InputComment,4);
      }
      elsif ($InputLine eq "XMAXOFFSET") {
        ($xMaxOffset[0],$xMaxOffset[1],$xMaxOffset[2],$InputComment)=split(' ',$InputComment,4);
      }
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName) (inputline=$InputLine)\n";
      }
    }
    
    ampsConfigLib::ChangeValueOfArray("double Orbiter::DomainSize::xMinOffset\\[3\\]",\@xMinOffset,"main/Orbiter.cpp");
    ampsConfigLib::ChangeValueOfArray("double Orbiter::DomainSize::xMaxOffset\\[3\\]",\@xMaxOffset,"main/Orbiter.cpp");
  } 
  elsif ($InputLine eq "CALCULATEDRAGCOEFFICIENT") {
    my $flag;
    
    ($flag,$InputLine)=split(' ',$InputComment,2);
    
    if ($flag eq "ON") {
      ampsConfigLib::ChangeValueOfVariable("bool Orbiter::Sampling::DragCoefficient::SamplingMode","true","main/Orbiter.cpp");
    }
    elsif ($flag eq "OFF") {
      ampsConfigLib::ChangeValueOfVariable("bool Orbiter::Sampling::DragCoefficient::SamplingMode","false","main/Orbiter.cpp");
    }
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }    
  } 
    
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
  
}
