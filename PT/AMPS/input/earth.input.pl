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
 
  $InputLine=~s/[:,;=]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  
  if ($InputLine eq "MESHSIGNATURE") {
    $line=~s/[=():]/ /g;
    ($InputLine,$line)=split(' ',$line,2);
    ($InputLine,$line)=split(' ',$line,2);

    ampsConfigLib::ChangeValueOfVariable("char Earth::Mesh::sign\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputLine."\"","main/Earth.cpp");
  }
  elsif ($InputLine eq "SEP") {
    ($s0,$s1)=split(' ',$InputComment,2);

    if ($s0 eq "ON") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_SEP__MODE_","_PIC_MODE_ON_","main/Earth.dfn");
    } 
    elsif($s0 eq "OFF") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_SEP__MODE_","_PIC_MODE_OFF_","main/Earth.dfn");
    }
    else {
      die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
    }    
  }
  elsif ($InputLine eq "GCR") {
    ($s0,$s1)=split(' ',$InputComment,2);

    if ($s0 eq "ON") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_GCR__MODE_","_PIC_MODE_ON_","main/Earth.dfn");
    } 
    elsif($s0 eq "OFF") {
      ampsConfigLib::RedefineMacro("_PIC_EARTH_GCR__MODE_","_PIC_MODE_OFF_","main/Earth.dfn");
    }
    else {
      die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
    }     
  }
  
  #sampling of the gyro-radius and gyro-frecuency in the computational domain
  elsif ($InputLine eq "PARTICLEDATASAMPLINGMODE") {
    ($s0,$InputComment)=split(' ',$InputComment,2);

    if ($s0 eq "ON") {
      ampsConfigLib::ChangeValueOfVariable("bool Earth::Sampling::ParticleData::SamplingMode","true","main/Earth_Sampling.cpp");      
      ($s0,$InputComment)=split(' ',$InputComment,2);
      
      if ($s0 eq "FUNCTION") {        
        ($s1,$s0)=split('!',$line,2);
        $s1=~s/[,;=]/ /g;
        
        ($s0,$s1)=split(' ',$s1,2);             
        ($s0,$s1)=split(' ',$s1,2);                
        ($s0,$s1)=split(' ',$s1,2);        
        ($s0,$s1)=split(' ',$s1,2);
                
        ampsConfigLib::AddLine2File("#ifdef _PIC_USER_DEFING_PARTICLE_SAMPLING__NODE_ \n#undef _PIC_USER_DEFING_PARTICLE_SAMPLING__NODE_ \n#endif \n\n#define _PIC_USER_DEFING_PARTICLE_SAMPLING__NODE_(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s,node)   ".$s0."(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s,node)\n","pic/picGlobal.dfn");        
      }
      else {
        die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }        
    }
    elsif ($s0 eq "OFF") {
      ampsConfigLib::ChangeValueOfVariable("bool Earth::Sampling::ParticleData::SamplingMode","false","main/Earth_Sampling.cpp");
    }
    else {
      die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
    }
  }
  
  #the number, locations, energy range, and the number of the ebergy intervals used in the spherical sampling surfaces 
  elsif ($InputLine eq "SPHERICALSHELLS") {
    ($s0,$InputComment)=split(' ',$InputComment,2);
    
    if ($s0 eq "ON")  {
      my @ShellRadiiTable;
      my $ShellRadiiTableLength=0;
      
      #sampling will occurs
      ampsConfigLib::ChangeValueOfVariable("bool Earth::Sampling::SamplingMode","true","main/Earth_Sampling.cpp");
      
      #read the rest of the line
      while (defined $InputComment) {
        ($s0,$InputComment)=split(' ',$InputComment,2);
      
        #cgeck wether the entry is the list of the shells radii
        if ($s0 eq "X") {
          while (defined $InputComment) {
            ($s0,$InputComment)=split(' ',$InputComment,2);
            
            if ( ($s0 eq "EMIN") || ($s0 eq "EMAX") || ($s0 eq "NLEVELS") ) {
              #the entry is a keyword of the section -> exist and parse the rest of the line
              last;
            }
            else {
              #the entry seems to be the next location of the samplign shell -> add it the list of the radii
              push(@ShellRadiiTable,$s0);
              $ShellRadiiTableLength++;
            }            
          }
          
          #the list of the shell's radii is complete -> add it to the sources 
          ampsConfigLib::ChangeValueOfArray("double Earth::Sampling::SampleSphereRadii\\[Earth::Sampling::nSphericalShells\\]",\@ShellRadiiTable,"main/Earth_Sampling.cpp");
          ampsConfigLib::ChangeValueOfVariable("    const int nSphericalShells",$ShellRadiiTableLength,"main/Earth.h");                   
        }
        
        #check whether the entry is another setting parameter
        if ($s0 eq "EMIN") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("double Earth::Sampling::Fluency::minSampledEnergy",$s0,"main/Earth.cpp");         
        }
        elsif ($s0 eq "EMAX") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("double Earth::Sampling::Fluency::maxSampledEnergy",$s0,"main/Earth.cpp");          
        }
        elsif ($s0 eq "NLEVELS") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("int Earth::Sampling::Fluency::nSampledLevels",$s0,"main/Earth.cpp");
        }
        else {
          die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
        }       
      }      
    }
    elsif ($s0 eq "OFF") {
       ampsConfigLib::ChangeValueOfVariable("bool Earth::Sampling::SamplingMode","false","main/Earth_Sampling.cpp");
    }     
    else {
      die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
    }  
  }
  
  #parameters of the impulse source 
  elsif ($InputLine eq "IMPULSESOURCE") {
    ($s0,$InputComment)=split(' ',$InputComment,2);
    
    if ($s0 eq "ON") {
      #the model is turned "ON" -> Set up the model parameters 
      my $nTotalInjectionLocations=0;
      my (@x0All,@x1All,@x2All,@SpecAll,@TimeAll,@SourceAll,@nPartAll);
            
      ampsConfigLib::ChangeValueOfVariable("bool Earth::ImpulseSource::Mode","true","main/ImpulseSource.cpp"); 
            
      while (defined $InputComment) {
        ($s0,$InputComment)=split(' ',$InputComment,2);
        
        if ($s0 eq "NEW") {
          $nTotalInjectionLocations++;
        }
        elsif ($s0 eq "X") {
          my ($x0,$x1,$x2);
          
          ($x0,$x1,$x2,$InputComment)=split(' ',$InputComment,4);
          
          push(@x0All,$x0);
          push(@x1All,$x1);
          push(@x2All,$x2);          
        }
        elsif ($s0 eq "SPEC") {
          my $nspec;
          
          ($s0,$InputComment)=split(' ',$InputComment,2);
          $nspec=ampsConfigLib::GetElementNumber($s0,\@SpeciesList);
          
          push(@SpecAll,$nspec);
        }
        elsif ($s0 eq "TIME") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          push(@TimeAll,$s0);
        }
        elsif ($s0 eq "NPART") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          push(@nPartAll,$s0);
        }        
        elsif ($s0 eq "SOURCE") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          push(@SourceAll,$s0);
        }
        
        elsif ($s0 eq "SPECTRUM") {
          $InputComment=~s/[()]/ /;
          $InputComment=~s/[()]/ /;
          
          ($s0,$InputComment)=split(' ',$InputComment,2);
          
          if ($s0 eq "CONSTANT") {
            ($s0,$InputComment)=split(' ',$InputComment,2);
            ampsConfigLib::ChangeValueOfVariable("double Earth::ImpulseSource::EnergySpectrum::Constant::e",$s0,"main/ImpulseSource.cpp"); 
            ampsConfigLib::ChangeValueOfVariable("int Earth::ImpulseSource::EnergySpectrum::Mode","Earth::ImpulseSource::EnergySpectrum::Mode_Constatant","main/ImpulseSource.cpp");           
          }
          else {
            die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
          }
          
        }
        
        else {
          die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
        }
      }
      
      #in case locations are defined -> set the model
      if ($nTotalInjectionLocations != 0) {
        #construct the list for the inpulse source data structure 
        my ($TotalDataStructureLine,$SingleEntry,$i);
        
        #{spec,time,false,{x0,x1,x2},Source};
        for ($i=0;$i<$nTotalInjectionLocations;$i++) {
          $SingleEntry="{".$SpecAll[$i].",".$TimeAll[$i].",false,{".$x0All[$i].",".$x1All[$i].",".$x2All[$i]."},".$SourceAll[$i].",".$nPartAll[$i]."}";
          
          if ($i==0) {
            $TotalDataStructureLine="{".$SingleEntry;
          }
          else {
            $TotalDataStructureLine=$TotalDataStructureLine.",".$SingleEntry;
          }          
        }
        
        $TotalDataStructureLine=$TotalDataStructureLine."}";
        ampsConfigLib::ChangeValueOfVariable("Earth::ImpulseSource::cImpulseSourceData Earth::ImpulseSource::ImpulseSourceData\\[\\]",$TotalDataStructureLine,"main/ImpulseSource.cpp");
        ampsConfigLib::ChangeValueOfVariable("int Earth::ImpulseSource::nTotalSourceLocations",$nTotalInjectionLocations,"main/ImpulseSource.cpp");
      }
       
    }  
  }
  
  
  
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
  
}
