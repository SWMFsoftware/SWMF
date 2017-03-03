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
  
  #point particle source 
  elsif ($InputLine eq "POINTSOURCE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "ON") {
      #read the point source section
      my $SourceLocationCounter=0;
      my (@spec,@xLocation,@yLocation,@zLocation,@Temperature,@SourceRate);
      
      $InputComment=~s/[:]/ /g;
      
      while (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "SOURCE") {
          $SourceLocationCounter++;
        }
        elsif ($InputLine eq "SPEC") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @spec,$InputLine;
        }
        elsif ($InputLine eq "RATE") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @SourceRate,$InputLine          
        }
        elsif ($InputLine eq "TEMPERATURE") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @Temperature,$InputLine          
        }
        elsif ($InputLine eq "X") {
          my ($x,$y,$z);
          
          ($x,$y,$z,$InputComment)=split(' ',$InputComment,4);
          push @xLocation,$x; 
          push @yLocation,$y;
          push @zLocation,$z;
        }
        else {
          die "Option is unknown in the 'Point Source' section, line=$InputFileLineNumber, option=$InputLine,  ($InputFileName)\n";
        }       
      }
      
      #insert the point source data into the code 
      my $InjectionDataTable="";
      
      for (my $i=0;$i<$SourceLocationCounter;$i++) {
        if ($i!=0) {
          $InjectionDataTable=$InjectionDataTable.",";
        }        
        
        $InjectionDataTable=$InjectionDataTable."{".$SourceRate[$i].",".$Temperature[$i].",{".$xLocation[$i].",".$yLocation[$i].",".$zLocation[$i]."},".$spec[$i]."}";        
      }
      
      ampsConfigLib::ChangeValueOfVariable("int Orbiter::InjectionModel::PointSource::InjectionDataTableLength",$SourceLocationCounter,"main/ParticleInjection.cpp"); 
      ampsConfigLib::ChangeValueOfVariable("Orbiter::InjectionModel::PointSource::cInjectionData Orbiter::InjectionModel::PointSource::InjectionDataTable\\[\\]","{".$InjectionDataTable."}","main/ParticleInjection.cpp");             
    }
    else {      
      ampsConfigLib::ChangeValueOfVariable("int Orbiter::InjectionModel::PointSource::InjectionDataTableLength","0","main/ParticleInjection.cpp"); 
    }
  }
  
  #source of the model particles due to the particle ejectino from faces 
  elsif ($InputLine eq "FACEEJECTION") {    
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "ON") {
      #read the point source section
      my $SourceLocationCounter=0;
      my (@spec,@faceat,@Temperature,@SourceRate);
      
      $InputComment=~s/[:]/ /g;
      
      while (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "SOURCE") {
          $SourceLocationCounter++;
        }
        elsif ($InputLine eq "SPEC") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @spec,$InputLine;
        }
        elsif ($InputLine eq "RATE") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @SourceRate,$InputLine          
        }
        elsif ($InputLine eq "TEMPERATURE") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @Temperature,$InputLine          
        }
        elsif ($InputLine eq "FACEAT") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @faceat,$InputLine           
        }
        else {
          die "Option is unknown in the 'Point Source' section, line=$InputFileLineNumber, option=$InputLine,  ($InputFileName)\n";
        }       
      }

      #insert the point source data into the code 
      my $InjectionDataTable="";
      
      for (my $i=0;$i<$SourceLocationCounter;$i++) {
        if ($i!=0) {
          $InjectionDataTable=$InjectionDataTable.",";
        }        
        
        $InjectionDataTable=$InjectionDataTable."{".$faceat[$i].",".$SourceRate[$i].",".$Temperature[$i].",".$spec[$i].",NULL}";        
      }
  
      ampsConfigLib::ChangeValueOfVariable("int Orbiter::InjectionModel::FaceEjection::InjectionDataTableLength",$SourceLocationCounter,"main/ParticleInjection.cpp"); 
      ampsConfigLib::ChangeValueOfVariable("Orbiter::InjectionModel::FaceEjection::cInjectionData Orbiter::InjectionModel::FaceEjection::InjectionDataTable\\[\\]","{".$InjectionDataTable."}","main/ParticleInjection.cpp");             
    }
    else {
      #the injection model is turned off 
      ampsConfigLib::ChangeValueOfVariable("int Orbiter::InjectionModel::FaceEjection::InjectionDataTableLength","0","main/ParticleInjection.cpp");
    }    
  }
  
  #parameters of the upstream gas injection source 
  elsif ($InputLine eq "UPSTREAMSOURCEMODE") {    
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "ON") {
      ampsConfigLib::ChangeValueOfVariable("bool Orbiter::UpstreamBC::UpstreamSourceMode","true","main/Orbiter.cpp"); 
    }
    elsif ($InputLine eq "OFF") {
      ampsConfigLib::ChangeValueOfVariable("bool Orbiter::UpstreamBC::UpstreamSourceMode","false","main/Orbiter.cpp");
    }
    else {
      die "Option is unknown in the 'Point Source' section, line=$InputFileLineNumber, option=$InputLine,  ($InputFileName)\n";
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
