#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
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


my $InputFileName=$ARGV[0]; # moon.input.Assembled.Block;
my $SpeciesFileName = $InputFileName; $SpeciesFileName =~ s/\.Block$/.Species/;

my $WorkingSourceDirectory=$ARGV[1];  # srcTemp
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
my $s3;

#source process id counter
my $SourceProcessID=0;

#the line symbolic id of the source processes
my @SourceProcessesSymbolicID;

#Modification of the surface content by the source processes
my @SourceModifySurfaceSpeciesAbundance;

#markers
my $MARKER__CALCULATE_SOURCE_FLUX_WITH_USER_DEFINED_FUNCTIONS; #substitude marker for calculation of the source rates for the user defined  source processes in the Exospehre.cpp % $MARKER:CALCULATE-SOURCE-FLUX-WITH-USER-DEFINED-FUNCTIONS$
my $MARKER__GENERATE_PARTICLE_PROPERTIES_WITH_USER_DEFINED_FUNCTIONS;#substitute marker $MARKER:GENERATE-PARTICLE-PROPERTIES-WITH-USER-DEFINED-FUNCTIONS$ for generation of the properties of the injected particles
my $MARKER__RESERVE_CELL_SAMPLING_DATA_BUFFER;#reserve space for sampling local density that is due to the user defined source process (Exosphere.cpp)
my $MARKER__USER_DEFINED_TOTAL_SOURCE_RATE;#calculate the total injection rate with the user defined source processes (Exosphere.cpp);




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
open (EXOSPHERE_USER_DEFINITIONS,">>$WorkingSourceDirectory/main/UserDefinition.Exosphere.h") || die "Cannot opne file $WorkingSourceDirectory/main/UserDefinition.Exosphere.h\n";

print EXOSPHERE_USER_DEFINITIONS "//The following changes are added while reading input file. [".strftime("%m/%d/%Y", localtime)."]\n";

while ($line=<InputFile>) {
  ($InputFileLineNumber,$FileName)=split(' ',$line);
  $line=<InputFile>;
  
  $LineOriginal=$line;
  chomp($line);

  for (my $n=1;$n<length($line)-1;$n++) {
    if (substr($line,$n,1) eq ":") {
      if ( (substr($line,$n-1,1) ne ":") && (substr($line,$n+1,1) ne ":")) {
        my $t=$line;
        $line=substr($line,0,$n)." ".substr($line,$n+1,length($line)-$n);
      }
    }
  }

  ($InputLine,$InputComment)=split('!',$line,2);
  $InputLine=uc($InputLine);
  chomp($InputLine);
 
  $InputLine=~s/=/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  
  #read the line 
  if ($InputLine eq "PHOTOLYTICREACTIONS") {
    my $ReactionFlag=0;
    my $ReactionProcessor;
    my $LifeTimeFunction;
    my $t;
    
    $line=~s/=/ /g;
    ($InputLine,$line)=split(' ',$line,2);
    
    while (defined $InputComment) {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ($t,$line)=split(' ',$line,2);
      $InputLine=~s/ //g;
      
      if ($InputLine eq "ON") {
        $ReactionFlag=1;
        ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS_MODE_\n#define _PIC_PHOTOLYTIC_REACTIONS_MODE_ _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_\n\n","main/UserDefinition.PIC.h");
      }
      elsif ($InputLine eq "OFF") {
        $ReactionFlag=0;
        ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS_MODE_\n#define _PIC_PHOTOLYTIC_REACTIONS_MODE_ _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_\n\n","main/UserDefinition.PIC.h");
      }
      elsif ($InputLine eq "REACTIONPROCESSOR") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        ($ReactionProcessor,$line)=split(' ',$line,2);
        $ReactionProcessor=~s/ //g;
        ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_\n#define _PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_(t0,t1,t2,t3,t4) $ReactionProcessor(t0,t1,t2,t3,t4);\n\n","main/UserDefinition.PIC.h");
      }
      elsif ($InputLine eq "LIFETIME") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        ($LifeTimeFunction,$line)=split(' ',$line,2);
        $LifeTimeFunction=~s/ //g;
        ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_\n#define _PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(t0,t1,t2,t3) $LifeTimeFunction(t0,t1,t2,t3);\n\n","main/UserDefinition.PIC.h");
      }
      else {
        die "Cannot recognize the option, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
  }
  
  
  elsif ($InputLine eq "SPICE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "ON") {
      ampsConfigLib::RedefineMacro("_EXOSPHERE__ORBIT_CALCUALTION__MODE_","_PIC_MODE_ON_","models/exosphere/Exosphere.h");
    }
    elsif ($InputLine eq "OFF") {
      ampsConfigLib::RedefineMacro("_EXOSPHERE__ORBIT_CALCUALTION__MODE_","_PIC_MODE_OFF_","models/exosphere/Exosphere.h");
    }
    else {
      die "Cannot recognize the option, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  
  elsif ($InputLine eq "TYPICALSOLARWINDCONDITIONS") {
    my @B=(0.0,0.0,0.0);
    my @v=(0.0,0.0,0.0);
    my $n=0.0;
    my $T=0.0;
    
    $InputComment=~s/[(),]/ /g;
     
    while (defined $InputComment) {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
      if ($InputLine eq "B") {
        ($B[0],$B[1],$B[2],$InputComment)=split(' ',$InputComment,4);
        $B[0]=~s/ //g;
        $B[1]=~s/ //g;
        $B[2]=~s/ //g;  
      }
      elsif ($InputLine eq "V") {
        ($v[0],$v[1],$v[2],$InputComment)=split(' ',$InputComment,4);
        $v[0]=~s/ //g;
        $v[1]=~s/ //g;
        $v[2]=~s/ //g;  
      }
      elsif ($InputLine eq "N") {
        ($n,$InputComment)=split(' ',$InputComment,2);
        $n=~s/ //g; 
      }      
      elsif ($InputLine eq "T") {
        ($T,$InputComment)=split(' ',$InputComment,2);
        $T=~s/ //g; 
      }
      else {
        die "The option is not found, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    
    #add the parameters to the source code
    ampsConfigLib::ChangeValueOfArray("static const double Exosphere_swVelocity_Typical\\[\\]",\@v,"models/exosphere/Exosphere.h");
    ampsConfigLib::ChangeValueOfArray("static const double Exosphere_swB_Typical\\[\\]",\@B,"models/exosphere/Exosphere.h");
    ampsConfigLib::ChangeValueOfVariable("static const double Exosphere_swTemperature_Typical",$T,"models/exosphere/Exosphere.h");
    ampsConfigLib::ChangeValueOfVariable("static const double Exosphere_swNumberDensity_Typical",$n,"models/exosphere/Exosphere.h");
    
  }
  elsif ($InputLine eq "ACCOMMODATIONCOEFFICIENT") {
    ($s0,$InputComment)=split(' ',$InputComment,2);
    
    if ($s0 eq "CONSTANT") {

      
      my @InputValue=(0)x$TotalSpeciesNumber;
      
      $InputComment=~s/[()=,]/ /g;

      while (defined $InputComment) {
        my $spec;
        my $value;
        
        ($s0,$InputComment)=split(' ',$InputComment,2);
        
        if ($s0 eq "CONST") {
          ($spec,$value,$InputComment)=split(' ',$InputComment,3);
          
          for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
            if ($SpeciesList[$i] eq $spec) {
              $InputValue[$i]=$value;
            }
          }
        }
        
      }
      
      #insert the array into the source code
      
      ampsConfigLib::ChangeValueOfArray("static const double AccomodationCoefficient\\[\\]",\@InputValue,"models/exosphere/Exosphere.h"); 
    }
    else {
      die "Error: not implemented, line=$InputFileLineNumber ($InputFileName)\n";
    }
    
  }
  
  elsif ($InputLine eq "SIMULATIONSTARTTIME") {
    chomp($LineOriginal);
    ($s0,$s1)=split('!',$LineOriginal,2);
    $s0=~s/=/ /g;
    ($s0,$s1)=split(' ',$s0,2);
    
    ampsConfigLib::ChangeValueOfVariable("static const char SimulationStartTimeString\\[\\]","\"$s1\"","models/exosphere/Exosphere.h");
  }
  
  elsif ($InputLine eq "ADDPHYSICALMODELHEADER") {  
    chomp($LineOriginal);
    $LineOriginal=~s/[(),=]/ /g;
    ($s0,$s1)=split('!',$LineOriginal,2);
    ($s0,$s1)=split(' ',$s0,2);

    while (defined $s1) {
      ($s0,$s1)=split(' ',$s1,2);
      ampsConfigLib::AddLine2File("#include \"$s0\"","main/UserDefinition.PIC.PhysicalModelHeaderList.h");
    }    
  }
  
  elsif ($InputLine eq "SURFACEVOLATILEDENSITY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "FUNCTION") {
      my $FunctionName;
      
      while (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "FUNCTION") {
          ($FunctionName,$InputComment)=split(' ',$InputComment,2);
          
          ampsConfigLib::AddLine2File("#undef _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_\n#define _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_(spec,el) $FunctionName\n\n","main/UserDefinition.Exosphere.h");
          ampsConfigLib::AddLine2File("#undef _EXOSPHERE__SURFACE_CONTENT_\n#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_\n\n","main/UserDefinition.Exosphere.h");
        }
      }
    }
    elsif ($InputLine eq "CONST") {
      my @SurfaceVolatileDensity=(0)x$TotalSpeciesNumber;
      my ($spec,$dens);
      
      $InputComment=~s/[()=,]/ /g;
      
      while (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "CONST") {
          ($spec,$dens,$InputComment)=split(' ',$InputComment,3);
          
          for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
            if ($SpeciesList[$i] eq $spec) {
              $SurfaceVolatileDensity[$i]=$dens;
            }
          }     
        }
      }
      
      
      my $OutputLine;
      
      $OutputLine="#ifndef _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED__DENSITY_TABLE_\n#define _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED__DENSITY_TABLE_\n";
      $OutputLine=$OutputLine."static const double ExosphereSurfaceContent_UserDefined_DensityTable[$TotalSpeciesNumber]={";
      
      for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
        $OutputLine=$OutputLine."$SurfaceVolatileDensity[$i]";
        
        if ($i ne $#SurfaceVolatileDensity) {
          $OutputLine=$OutputLine.",";
        }       
      }
      
      $OutputLine=$OutputLine."};\n";     
      $OutputLine=$OutputLine."#undef _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_\n#define _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_(spec,el) (ExosphereSurfaceContent_UserDefined_DensityTable[spec])\n";
      $OutputLine=$OutputLine."#undef _EXOSPHERE__SURFACE_CONTENT_\n#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_\n";
      $OutputLine=$OutputLine."#endif\n\n";
      
      ampsConfigLib::AddLine2File($OutputLine,"main/UserDefinition.Exosphere.h");
    }    
    elsif ($InputLine eq "FLUXBALANCE") {
      ampsConfigLib::AddLine2File("#undef _EXOSPHERE__SURFACE_CONTENT_\n#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__BALANCE_FLUXES_\n\n","main/UserDefinition.Exosphere.h");
    }    
    elsif ($InputLine eq "UNIFORM") {
      ampsConfigLib::AddLine2File("#undef _EXOSPHERE__SURFACE_CONTENT_\n#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__UNIFORM_\n\n","main/UserDefinition.Exosphere.h");
    }    
    elsif ($InputLine eq "RADIALDISTRIBUTION") {
      ampsConfigLib::AddLine2File("#undef _EXOSPHERE__SURFACE_CONTENT_\n#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__RADIAL_DISTRIBUTION_\n\n","main/UserDefinition.Exosphere.h");
    }       
    else {
      die "Cannot recognize the option, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }  
  elsif ($InputLine eq "SOURCE") {
    $InputComment=~s/[=(),]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    
    ($s0,$InputComment)=split(' ',$InputComment,2);
    $s0=~s/ //g;
    
    next if ($s0 ne "ON");
    
    
    if ($InputLine eq "PHOTONSTIMULATEDDESORPTION") {
      my $PhotonFlux_1AU=0;
      my @CrossSection=(0)x$TotalSpeciesNumber;
      my @minInjectionEnergy=(0)x$TotalSpeciesNumber;
      my @maxInjectionEnergy=(0)x$TotalSpeciesNumber;
      my $nspec;
      
      while (defined $InputComment) {
        ($InputLine,$s0,$InputComment)=split(' ',$InputComment,3);
        $InputLine=~s/ //g;
        $s0=~s/ //g;
        
        if ($InputLine eq "PHOTONFLUX_1AU") {
          $PhotonFlux_1AU=$s0;
        }
        elsif ($InputLine eq "CROSSSECTION") {
          ($s1,$InputComment)=split(' ',$InputComment,2);
          $s1=~s/ //g;
          $nspec=getSpeciesNumber($s0);
          $CrossSection[$nspec]=$s1;
        }
        elsif ($InputLine eq "INJECTIONVELOCITYRANGE") {
          $nspec=getSpeciesNumber($s0);
          ($s0,$s1,$InputComment)=split(' ',$InputComment,3);
          $s0=~s/ //g;
          $s1=~s/ //g;
          
          $minInjectionEnergy[$nspec]="_".$SpeciesList[$nspec]."__MASS_*pow($s0,2)/2.0";
          $maxInjectionEnergy[$nspec]="_".$SpeciesList[$nspec]."__MASS_*pow($s1,2)/2.0";
        }
        else {
          die "Cannot recognize the option, line=$InputFileLineNumber ($InputFileName)\n";
        }
      }
      
      #add the parameters to the source code
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_","_EXOSPHERE_SOURCE__ON_","models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfVariable("static const double PhotonStimulatedDesorption_PhotonFlux_1AU",$PhotonFlux_1AU,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double PhotonStimulatedDesorption_CrossSection\\[\\]",\@CrossSection,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double PhotonStimulatedDesorption_minInjectionEnergy\\[\\]",\@minInjectionEnergy,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double PhotonStimulatedDesorption_maxInjectionEnergy\\[\\]",\@maxInjectionEnergy,"models/exosphere/Exosphere.h");
      
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_",$SourceProcessID,"models/exosphere/Exosphere.h");
 
      push(@SourceModifySurfaceSpeciesAbundance,'true');
      push(@SourceProcessesSymbolicID,"\"PhotonStimulatedDesorption\"");
      $SourceProcessID++;
      
    }
    elsif ($InputLine eq "THERMALDESORPTION") {
      my @uThermal=(0)x$TotalSpeciesNumber;
      my @VibrationalFrequency=(0)x$TotalSpeciesNumber;
      my $nspec;
      
      while (defined $InputComment) {
        ($InputLine,$s0,$s1,$InputComment)=split(' ',$InputComment,4);
        $InputLine=~s/ //g;
        $s0=~s/ //g;
        $s1=~s/ //g;
        
        $nspec=getSpeciesNumber($s0);
        
        if ($InputLine eq "UTHERMAL") {
          $uThermal[$nspec]=$s1;
        }
        elsif ($InputLine eq "VIBRATIONALFREQUENCY") {
          $VibrationalFrequency[$nspec]=$s1;
        }
        else {
          die "Cannot recognize the option, line=$InputFileLineNumber ($InputFileName)\n";
        }
      }
      
      #add the parameters to the source code
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__THERMAL_DESORPTION_","_EXOSPHERE_SOURCE__ON_","models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double ThermalDesorption_uThermal\\[\\]",\@uThermal,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double ThermalDesorption_VibrationalFrequency\\[\\]",\@VibrationalFrequency,"models/exosphere/Exosphere.h");
      
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_",$SourceProcessID,"models/exosphere/Exosphere.h");     
      
      push(@SourceModifySurfaceSpeciesAbundance,'true');
      push(@SourceProcessesSymbolicID,"\"ThermalDesorption\"");
      $SourceProcessID++;
    }
    elsif ($InputLine eq "SPUTTERING") {
      my @Yield=(0)x$TotalSpeciesNumber;
      my @minInjectionEnergy=(0)x$TotalSpeciesNumber;
      my @maxInjectionEnergy=(0)x$TotalSpeciesNumber;
      
      while (defined $InputComment) {
        ($InputLine,$s0,$InputComment)=split(' ',$InputComment,3);
        $InputLine=~s/ //g;
        $s0=~s/ //g;
        
        if ($InputLine eq "YIELD") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          $InputLine=~s/ //g;
          
          $Yield[getSpeciesNumber($s0)]=$InputLine;
        }
        elsif ($InputLine eq "INJECTIONVELOCITYRANGE") {
          my $nspec=getSpeciesNumber($s0);
          ($s0,$s1,$InputComment)=split(' ',$InputComment,3);
          $s0=~s/ //g;
          $s1=~s/ //g;
          
          $minInjectionEnergy[$nspec]="_".$SpeciesList[$nspec]."__MASS_*pow($s0,2)/2.0";
          $maxInjectionEnergy[$nspec]="_".$SpeciesList[$nspec]."__MASS_*pow($s1,2)/2.0";
        }
        else {
          die "Cannot recognize the option, line=$InputFileLineNumber ($InputFileName)\n";
        }        
      }
      
      #add the parameters to the source code
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_","_EXOSPHERE_SOURCE__ON_","models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double SolarWindSputtering_Yield\\[\\]",\@Yield,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double SolarWindSputtering_minInjectionEnergy\\[\\]",\@minInjectionEnergy,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double SolarWindSputtering_maxInjectionEnergy\\[\\]",\@maxInjectionEnergy,"models/exosphere/Exosphere.h");
      
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_",$SourceProcessID,"models/exosphere/Exosphere.h");
      
      push(@SourceModifySurfaceSpeciesAbundance,'true');
      push(@SourceProcessesSymbolicID,"\"SolarWindSputtering\"");
      $SourceProcessID++;
    }
    elsif ($InputLine eq "IMPACTVAPORIZATION") {
      my $HeliocentricDistance="_AU_";
      my $SourceRatePowerIndex=0;
      my @SourceRate=(0)x$TotalSpeciesNumber;
      my @SourceTemperature=(0)x$TotalSpeciesNumber;
      
      while (defined $InputComment) {
        ($InputLine,$s0,$InputComment)=split(' ',$InputComment,3);
        $InputLine=~s/ //g;
        $s0=~s/ //g;
        
        if ($InputLine eq "HELIOCENTRICDISTANCE") {
          $HeliocentricDistance=$s0;
        }
        elsif ($InputLine eq "SOURCERATEPOWERINDEX") {
          $SourceRatePowerIndex=$s0;
        }
        elsif ($InputLine eq "SOURCERATE") {
          ($s1,$InputComment)=split(' ',$InputComment,2);
          $s1=~s/ //g;
          $SourceRate[getSpeciesNumber($s0)]=$s1;
        }
        elsif ($InputLine eq "SOURCETEMPERATURE") {
          ($s1,$InputComment)=split(' ',$InputComment,2);
          $s1=~s/ //g;
          $SourceTemperature[getSpeciesNumber($s0)]=$s1;
        }
        else {
          die "Cannot recognize the option, line=$InputFileLineNumber ($InputFileName)\n";
        }
      }
      
      #add the parameters of the input file to the code      
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_","_EXOSPHERE_SOURCE__ON_","models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfVariable("static const double ImpactVaporization_HeliocentricDistance",$HeliocentricDistance,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfVariable("static const double ImpactVaporization_SourceRatePowerIndex",$SourceRatePowerIndex,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double ImpactVaporization_SourceRate\\[\\]",\@SourceRate,"models/exosphere/Exosphere.h");
      ampsConfigLib::ChangeValueOfArray("static const double ImpactVaporization_SourceTemeprature\\[\\]",\@SourceTemperature,"models/exosphere/Exosphere.h");

      ampsConfigLib::RedefineMacro("_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_",$SourceProcessID,"models/exosphere/Exosphere.h");
            
      push(@SourceModifySurfaceSpeciesAbundance,'false');
      push(@SourceProcessesSymbolicID,"\"ImpactVaposization\"");
      $SourceProcessID++;
    }
    elsif ($InputLine eq "USERDEFINED") {
      my $Code;
      my $SourceRate;
      my $GenerateParticleProperties;
      my $InitSurfaceSourceDistribution="";
      
      my $UserDefinedModifySurfaceAbundance;
      
      
      if ($s0 eq "ON") {
        #the user defined source is in use
        #parse the original line from the input file
        my $s0original;
        my $s1original;
        my $s2original;
        
        
        ($line,$s0original)=split('!',$line,2);
        
        $line=~s/=/ /g;
        ($s0,$s1,$s2,$line)=split(' ',$line,4);
        
        
        while (defined $InputComment) {
          ($s0,$s1,$InputComment)=split(' ',$InputComment,3);
          
          if (defined $InputComment) {
            $InputComment=~s/^\s+//; #remove spaces from the begining of the line
          }
          
          ($s0original,$s1original,$line)=split(' ',$line,3);
          
          if (defined $line) {
            $line=~s/^\s+//; #remove spaces from the begining of the line
          }
          
          if (!defined $s0) {
            $s0="nothing";
          }
          
          if ($s0 eq "SOURCEPROCESSCODE") {
            $Code=$s1original;
          }
          elsif ($s0 eq "SOURCERATE") {
            $SourceRate=$s1original;
          }
          elsif ($s0 eq "GENERATEPARTICLEPROPERTIES") {
            $GenerateParticleProperties=$s1original;
          }
          elsif ($s0 eq "INITSURFACESOURCEDISTRIBUTION") {
            $InitSurfaceSourceDistribution=$s1original;
          }
          elsif ($s0 eq "MODIFYSURFACESPECIESABUNDANCE") {
            $UserDefinedModifySurfaceAbundance="0";
            
            if ($s1 eq "TRUE") {
              $s1='true';
            }
            elsif ($s1 eq "FALSE") {
              $s1='false';
            }
            else {
              die "Cannot recognize option [file=$InputFileName, line=$InputFileLineNumber]\n";
            }
            
            push(@SourceModifySurfaceSpeciesAbundance,$s1);
          }
          else {
            die "Cannot recognize option [file=$InputFileName, line=$InputFileLineNumber]\n";
          }
        }
        
        if (!defined $UserDefinedModifySurfaceAbundance) {
          print "ModifySurfaceSpeciesAbundance is not defined [file=$InputFileName, line=$InputFileLineNumber]\n";
        }
        elsif (!defined $Code) {
          die "SourceProcessCode is not defined [file=$InputFileName, line=$InputFileLineNumber]\n";
        }
        elsif (!defined $SourceRate) {
          die "SourceRate is not defined [file=$InputFileName, line=$InputFileLineNumber]\n";
        }
        elsif (!defined $GenerateParticleProperties) {
          die "GenerateParticleProperties is not defined [file=$InputFileName, line=$InputFileLineNumber]\n";
        }
        
        
        
        my $SourceCode="_EXOSPHERE_SOURCE__ID__USER_DEFINED__".$SourceProcessID."_".$Code."_";
        
        print EXOSPHERE_USER_DEFINITIONS "#define $SourceCode  $SourceProcessID\n";
        print EXOSPHERE_USER_DEFINITIONS "#define _EXOSPEHRE_SOURCE__USER_DEFINED__".$SourceProcessID."_".$Code."_  _EXOSPHERE_SOURCE__ON_\n";
        print EXOSPHERE_USER_DEFINITIONS "\n#undef _EXOSPHERE__USER_DEFINED_SOURCE_MODEL__MODE_\n#define _EXOSPHERE__USER_DEFINED_SOURCE_MODEL__MODE_  _EXOSPHERE_SOURCE__ON_\n";
        
        $MARKER__CALCULATE_SOURCE_FLUX_WITH_USER_DEFINED_FUNCTIONS=$MARKER__CALCULATE_SOURCE_FLUX_WITH_USER_DEFINED_FUNCTIONS."\nFluxSourceProcess[$SourceCode]=$SourceRate(spec,SphereDataPointer);\n";
        
        if ($InitSurfaceSourceDistribution ne "") {
          $MARKER__CALCULATE_SOURCE_FLUX_WITH_USER_DEFINED_FUNCTIONS=$MARKER__CALCULATE_SOURCE_FLUX_WITH_USER_DEFINED_FUNCTIONS."if (FluxSourceProcess[$SourceCode]>0.0) $InitSurfaceSourceDistribution();\n";
        }
        
#       $MARKER__GENERATE_PARTICLE_PROPERTIES_WITH_USER_DEFINED_FUNCTIONS=$MARKER__GENERATE_PARTICLE_PROPERTIES_WITH_USER_DEFINED_FUNCTIONS."\nelse if (SourceProcessID==$SourceCode) {\nflag=$GenerateParticleProperties(spec,(PIC::ParticleBuffer::byte*)tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere);\nSourceProcessID=$SourceCode;\nif (flag==true) Sampling::CalculatedSourceRate[spec][$SourceCode]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;\n}\n";

        $MARKER__GENERATE_PARTICLE_PROPERTIES_WITH_USER_DEFINED_FUNCTIONS=$MARKER__GENERATE_PARTICLE_PROPERTIES_WITH_USER_DEFINED_FUNCTIONS."\nelse if (SourceProcessID==$SourceCode) {\nflag=$GenerateParticleProperties(spec,(PIC::ParticleBuffer::byte*)tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere);\n}\n"; 
        $MARKER__RESERVE_CELL_SAMPLING_DATA_BUFFER=$MARKER__RESERVE_CELL_SAMPLING_DATA_BUFFER."\nSamplingDensityOffset[$SourceCode]=CellSamplingDataOffset+SamplingLength;\nSamplingLength+=sizeof(double)*PIC::nTotalSpecies;\n";
        $MARKER__USER_DEFINED_TOTAL_SOURCE_RATE=$MARKER__USER_DEFINED_TOTAL_SOURCE_RATE."\nres+=$SourceRate(spec,SphereDataPointer);\n";
               
#        print "$SourceRate, $GenerateParticleProperties\n";
        
        push(@SourceProcessesSymbolicID,"\"$Code\"");
        $SourceProcessID++;
      }
    }
    
    
  }
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }
  else {
    chomp($LineOriginal);
    print "The option is unknown: (line=$InputFileLineNumber, file=$FileName)\n\"$LineOriginal\"\n";
    die;
  }
  
  
}


#the total number of sources
my $TotalNumberOfExosphericSources=$SourceProcessID;
my $MaxSourceIdNumber=$SourceProcessID;

if ($MaxSourceIdNumber ne 0) {
  $MaxSourceIdNumber-=1;
}


ampsConfigLib::RedefineMacro("_EXOSPHERE__SOURCE_TOTAL_NUMBER_",$TotalNumberOfExosphericSources,"models/exosphere/Exosphere.h");
ampsConfigLib::RedefineMacro("_EXOSPHERE__SOURCE_MAX_ID_VALUE_",$MaxSourceIdNumber,"models/exosphere/Exosphere.h");
ampsConfigLib::ChangeValueOfArray("static const char _EXOSPHERE__SOURCE_SYMBOLIC_ID_[][100]",\@SourceProcessesSymbolicID,"models/exosphere/Exosphere.h");

ampsConfigLib::ChangeValueOfArray("static const bool Source_DeplitSurfaceSpeciesAbundance_Flag\\[\\]",\@SourceModifySurfaceSpeciesAbundance,"models/exosphere/Exosphere.h");


#close the file with the user definitions 
close(EXOSPHERE_USER_DEFINITIONS);

#substitue markers in 'Exosphere.cpp'
my @FileContent;

if (defined $MARKER__CALCULATE_SOURCE_FLUX_WITH_USER_DEFINED_FUNCTIONS) {
  my $line;
  
  open (FILEIN,"<$WorkingSourceDirectory/models/exosphere/Exosphere.cpp") || die "Cannot open file $WorkingSourceDirectory/models/exosphere/Exosphere.cpp\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$WorkingSourceDirectory/models/exosphere/Exosphere.cpp") || die "Cannot open file $WorkingSourceDirectory/models/exosphere/Exosphere.cpp\n";
  
  foreach (@FileContent) {   
    $_=~s/\$MARKER:CALCULATE-SOURCE-FLUX-WITH-USER-DEFINED-FUNCTIONS\$/$MARKER__CALCULATE_SOURCE_FLUX_WITH_USER_DEFINED_FUNCTIONS/;
    $_=~s/\$MARKER:GENERATE-PARTICLE-PROPERTIES-WITH-USER-DEFINED-FUNCTIONS\$/$MARKER__GENERATE_PARTICLE_PROPERTIES_WITH_USER_DEFINED_FUNCTIONS/;
    $_=~s/\$MARKER:RESERVE-CELL-SAMPLING-DATA-BUFFER\$/$MARKER__RESERVE_CELL_SAMPLING_DATA_BUFFER/;
    $_=~s/\$MARKER:USER-DEFINED-TOTAL-SOURCE-RATE\$/$MARKER__USER_DEFINED_TOTAL_SOURCE_RATE/;
      
    print FILEOUT "$_";
  }
  
  
  close (FILEOUT);
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

=comment
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
#      print "\#define $Macro $Value\n";
      $_="\#define $Macro $Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

#======================================== Add a line to a file ===========================================
sub AddLine2File {
  my $fname=$_[1];
  my $newline=$_[0];
  
  open (EDITFILE,">>$WorkingSourceDirectory/$fname") || die "Cannot open file $WorkingSourceDirectory/$fname\n";
  print EDITFILE "$newline\n";
  close(EDITFILE);
}

=cut





