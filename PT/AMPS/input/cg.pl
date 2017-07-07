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


my $InputFileName=$ARGV[0]; #"cg.input.Assembled.Block";   #$ARGV[0];
my $SpeciesFileName=$InputFileName; $SpeciesFileName =~ s/\.Block$/.Species/; #"cg.input.Assembled.Species";
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
  $line=$InputLine;
  
  $InputLine=uc($InputLine);
  chomp($InputLine);
 
  $InputLine=~s/[=():,]/ /g;
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

  #change the starting locations at which tracking of the partucle trajecoties begis
  elsif ($InputLine eq "TRACINGSURFACERADIUS") { #TracingSurfaceRadius
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
     ampsConfigLib::ChangeValueOfVariable("const double TracingSurfaceRadius",$InputLine,"main/Comet.h"); 
  }
  
  #the range of the dust grain speed (used for evaluation of the local time step)
  elsif ($InputLine eq "DUSTVELOCITYLIMIT") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "MIN") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity",$InputLine,"main/Comet.cpp");   
    }
    elsif ($InputLine eq "MAX") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity",$InputLine,"main/Comet.cpp");   
    }     
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }     
  }
  
  #turn on sampling of the Rosina data
  elsif ($InputLine eq "ROSINADATASAMPLING") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if  ($InputLine eq "ON") {
      ampsConfigLib::RedefineMacro("_COMET_SAMPLE_ROSINA_DATA_","_PIC_MODE_ON_","main/Comet.dfn");   
    }
    elsif ($InputLine eq "OFF") {
      ampsConfigLib::RedefineMacro("_COMET_SAMPLE_ROSINA_DATA_","_PIC_MODE_OFF_","main/Comet.dfn");         
    }     
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }

  #define settingsof the End-of-Mission model
  elsif ($InputLine eq "ENDOFMISSION") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);

    ### EndOfMission::Test ####
    if ($InputLine eq "TEST") {
      #settings of the spherical test 
      ($InputLine,$InputComment)=split(' ',$InputComment,2);

      if ($InputLine eq "MODE") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);

        if ($InputLine eq "POINT") {
          ampsConfigLib::ChangeValueOfVariable("int SphericalNucleusTest_Mode","SphericalNucleusTestMode_Point","main/RosinaMeasurements_Liouville.cpp"); 
        }
        elsif ($InputLine eq "RANDOM") {
          ampsConfigLib::ChangeValueOfVariable("int SphericalNucleusTest_Mode","SphericalNucleusTestMode_Random","main/RosinaMeasurements_Liouville.cpp"); 
        }
        elsif ($InputLine eq "OFF") {
          ampsConfigLib::ChangeValueOfVariable("int SphericalNucleusTest_Mode","SphericalNucleusTestMode_OFF","main/RosinaMeasurements_Liouville.cpp"); 
        }
        else {
          warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
          die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
        }
      }
      elsif ($InputLine eq "SOURCERATE") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        ampsConfigLib::ChangeValueOfVariable("static double SphericalNucleusTest_SourceRate",$InputLine,"main/RosinaMeasurements_Liouville.cpp"); 
      }         
      elsif ($InputLine eq "TEMPERATURE") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        ampsConfigLib::ChangeValueOfVariable("static double SphericalNucleusTest_Temperature",$InputLine,"main/RosinaMeasurements_Liouville.cpp");
      }
      elsif ($InputLine eq "LOCATION") {
        my (@x,$i);
        
        for ($i=0;$i<3;$i++) {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @x, $InputLine;
        }
        
        ampsConfigLib::ChangeValueOfArray("static double SphericalNucleusTest_Location\\[3\\]",\@x,"main/RosinaMeasurements_Liouville.cpp");
      }
      elsif ($InputLine eq "LINEOFSIGHTRG") {
        my (@x,$i);
        
        for ($i=0;$i<3;$i++) {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @x, $InputLine;
        }
        
        ampsConfigLib::ChangeValueOfArray("static double SphericalNucleusTest_LineOfSightRG\\[3\\]",\@x,"main/RosinaMeasurements_Liouville.cpp");
      }  
      elsif ($InputLine eq "LINEOFSIGHTNG") {
        my (@x,$i);
        
        for ($i=0;$i<3;$i++) {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          push @x, $InputLine;
        }
        
        ampsConfigLib::ChangeValueOfArray("static double SphericalNucleusTest_LineOfSightNG\\[3\\]",\@x,"main/RosinaMeasurements_Liouville.cpp");
      }          
      else {
        warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
        die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    
    ### EndOfMission::SelfShadowingNG
    elsif ($InputLine eq "SELFSHADOWINGNG") {          
      while (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
        if ($InputLine eq "MODE") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          
          if ($InputLine eq "OFF") {
            ampsConfigLib::ChangeValueOfVariable("static const int SelfShadowingNGMode","SelfShadowingNG_ModeOff","main/RosinaMeasurements_Liouville.cpp");          
          }
          elsif ($InputLine eq "SIMPLE") {
            ampsConfigLib::ChangeValueOfVariable("static const int SelfShadowingNGMode","SelfShadowingNG_ModeSimple","main/RosinaMeasurements_Liouville.cpp");          
          }
          else {
            warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
            die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";          
          }
        }
        elsif ($InputLine eq "COSANGLELIMIT") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("static const double SelfShadowingNGCosAngleLimit",$InputLine,"main/RosinaMeasurements_Liouville.cpp");                  
        }
        else {
          warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
          die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";          
        }
      }
    }
    
    #### EndOfMission::Correction ###
    elsif ($InputLine eq "CORRECTION") {
      #settings of the correction model
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
            
      if ($InputLine eq "ADJUSTSURFACEINJECTIONRATE") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "ON") {
          ampsConfigLib::ChangeValueOfVariable("static bool AdjustSurfaceInjectionRate","true","main/RosinaMeasurements_Liouville.cpp");
        }
        else{
          ampsConfigLib::ChangeValueOfVariable("static bool AdjustSurfaceInjectionRate","false","main/RosinaMeasurements_Liouville.cpp");
        }
      }
      ### EndOfMission::Correction::NG ###
      elsif ($InputLine eq "NG") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
                
        if ($InputLine eq "MODE") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          
          if ($InputLine eq "OFF") {
            ampsConfigLib::ChangeValueOfVariable("static double NudeGaugeDensitySinCorrectionFactor","0.0","main/RosinaMeasurements_Liouville.cpp");
          }
          elsif ($InputLine eq "DENSITYSINCORRECTIONFACTOR") {
            ($InputLine,$InputComment)=split(' ',$InputComment,2);
            ampsConfigLib::ChangeValueOfVariable("static double NudeGaugeDensitySinCorrectionFactor",$InputLine,"main/RosinaMeasurements_Liouville.cpp");            
          }
          else {
            warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
            die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
          }
        }               
      }
      else {
        warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");        
        die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    
    ### EndOfMission::Limit ####
    elsif ($InputLine eq "LIMIT" ) {
      #settings of the limits 
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "NG") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "SOLIDANGLECUTOFF") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("static double CutoffNudeGaugeSolidAngle",$InputLine,"main/RosinaMeasurements_Liouville.cpp");                  
        }
        else {
          warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");          
          die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
        }
      }
      else {
        warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
        die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
      }

    }
    
    ### REFERENCE FILE CORRECTION
    elsif ($InputLine eq "REFERENCEFILECORRECTION" ) {
      my ($refname,$mode,$t);
      
      $mode="false";
      $line=~s/[=():,]/ /g;
      ($InputLine,$line)=split(' ',$line,2);
      ($InputLine,$line)=split(' ',$line,2);
      
      print "$line\n";
      
      while (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        ($t,$line)=split(' ',$line,2);
        
        if ($InputLine eq "MODE") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          ($t,$line)=split(' ',$line,2);
          
          if ($InputLine eq "ON") {
            $mode="true";
          }  
          else {
            $mode="false";
          }        
         
          ampsConfigLib::ChangeValueOfVariable("bool ApplyRefCorrection_$refname",$mode,"main/RosinaMeasurements_Liouville.cpp");
        }
        elsif ($InputLine eq "REF") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          ($refname,$line)=split(' ',$line,2);     
        }
        else {
          warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
          die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
        }    
      }            
    }
    
    ### UNSERTANTIES DUE TO LOCATION AND ORIENTATION OF THE SPACECRAFT ###
    elsif ($InputLine eq "TRAJECTORYUNSERTANTY") {
      my ($Mode,$nTests,$Radius,$AngularLimit);

      $Mode="false";
      $nTests=10;
      $Radius=50.0;
      $AngularLimit=5.0;

      while (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "NTESTS") { ($nTests,$InputComment)=split(' ',$InputComment,2);}
        elsif ($InputLine eq "RADIUS") {($Radius,$InputComment)=split(' ',$InputComment,2);}
        elsif ($InputLine eq "ANGULARLIMIT") {($AngularLimit,$InputComment)=split(' ',$InputComment,2);}
        elsif ($InputLine eq "MODE") {
          ($Mode,$InputComment)=split(' ',$InputComment,2);
          
          if ($Mode eq "ON") {$Mode="true";}
          elsif ($Mode eq "OFF") {$Mode="false";}
          else {
            warn("Option is unknown ($Mode), line=$InputFileLineNumber ($InputFileName)");
            die "Option is unknown ($Mode), line=$InputFileLineNumber ($InputFileName)\n";
          }
        }
      }

      #save trajectory unsertanty quantification parameters
      ampsConfigLib::ChangeValueOfVariable("int TrajectoryUncertanties_nTest",$nTests,"main/RosinaMeasurements_Liouville.cpp");
      ampsConfigLib::ChangeValueOfVariable("double TrajectoryUncertanties_Radius",$Radius,"main/RosinaMeasurements_Liouville.cpp");
      ampsConfigLib::ChangeValueOfVariable("bool TrajectoryUncertanties_Search",$Mode,"main/RosinaMeasurements_Liouville.cpp");
      ampsConfigLib::ChangeValueOfVariable("double TrajectoryUncertanties_AngularLimit",$AngularLimit,"main/RosinaMeasurements_Liouville.cpp");
    }
          
    ###INDIVIDUAL COMMANDS ###
    elsif ($InputLine eq "SIMULATIONDATAPOINTSTEP") {  
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("static int RosinaDataSimulationStep",$InputLine,"main/RosinaMeasurements_Liouville.cpp"); 
    }
    
    #disregard instrument orientation flag
    elsif ($InputLine eq "DISREGARDINSTRUMENTORIENTATION") {  
      ($InputLine,$InputComment)=split(' ',$InputComment,2);

      if ($InputLine eq "ON") {
        ampsConfigLib::ChangeValueOfVariable("bool DisregardInstrumentOrientationFlag","true","main/RosinaMeasurements_Liouville.cpp");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::ChangeValueOfVariable("bool DisregardInstrumentOrientationFlag","false","main/RosinaMeasurements_Liouville.cpp");
      }
      else {
        warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
        die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
      }
    }   

    #the number of the test in calcualtion of the solid angle occupied by the nucleus
    elsif ($InputLine eq "SOLIDANGLETEST") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::RedefineMacro("_ROSINA_SAMPLE__LIOUVILLE__GET_SOLID_ANGLE__TEST_NUMBER_","$InputLine","main/RosinaMeasurements.dfn");
    }

    #the number of the test points on each surface element when evaluating flux from the element that reaches the instrument
    elsif ($InputLine eq "LOCATIONTEST") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::RedefineMacro("_ROSINA_SAMPLE__LIOUVILLE__EVALUATE_LOCATION__TEST_NUMBER_","$InputLine","main/RosinaMeasurements.dfn");
    }
        
    elsif ($InputLine eq "INJECTIONDISTRIBUTIONMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "FLUX") {
        ampsConfigLib::ChangeValueOfVariable("static int VelocityInjectionMode","_ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VERTICAL_FLUX_MAXWELLIAN_","main/RosinaMeasurements_Liouville.cpp");                         
      }
      elsif ($InputLine eq "RANDOMFLUX") {
        ampsConfigLib::ChangeValueOfVariable("static int VelocityInjectionMode","_ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__RANDOMLY_DIRECTED_FLUX_MAXWELLAIN_","main/RosinaMeasurements_Liouville.cpp");
      }
      elsif ($InputLine eq "VELOCITY") {
        ampsConfigLib::ChangeValueOfVariable("static int VelocityInjectionMode","_ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VELOCITY_MAXWELLIAN_","main/RosinaMeasurements_Liouville.cpp");
      }      
      else {
        warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
        die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
      }
    } 
    else {
      warn("Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)");
      die "Option is unknown ($InputLine), line=$InputFileLineNumber ($InputFileName)\n";
    }
  
    ### END Of EndOfMission###
  }

  #redefine a value of a macro 
  elsif ($InputLine eq "DEFINE") {
    my ($macro,$value,$s0,$s1);
    
    ($InputLine,$InputComment)=split('!',$line,2);
    ($s0,$macro,$value,$s1)=split(' ',$InputLine,4);
    
    $s0=$macro;    
    $s0=~s/[()=]/ /g;
    ($s0,$s1)=split(' ',$s0,2);

    if (!defined $value) {
      $value=" ";
    }

    ampsConfigLib::AddLine2File("\n#undef $s0\n#define $macro $value\n","main/Comet.dfn");    
  }
  
  #define the mesh sigrature that will be used in the simuation (for generation of the particular mesh settings and loadging particular model data files)
  elsif ($InputLine eq "MESHSIGNATURE") {
    $line=~s/[=():]/ /g;
    ($InputLine,$line)=split(' ',$line,2);
    ($InputLine,$line)=split(' ',$line,2);
    
    ampsConfigLib::ChangeValueOfVariable("char Comet::Mesh::sign\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputLine."\"","main/Comet.cpp");
  }
  
  #mean density of the dust grains (DustMeanDensity)
  elsif ($InputLine eq "DUSTMEANDENSITY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2); 

    ampsConfigLib::ChangeValueOfVariable("double ElectricallyChargedDust::MeanDustDensity",$InputLine,"models/dust/Dust.cpp");
  } 

  #set the drag coefficient (GrainDragCoefficient)
  elsif ($InputLine eq "DUSTDRAGCOEFFICIENT") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);

    ampsConfigLib::ChangeValueOfVariable("double GrainDragCoefficient",$InputLine,"main/Comet.h");
  } 


  #type of the gas velocity distribution at the nucleus 
  elsif ($InputLine eq "INITIALGASVELOCITYDISTRIBUTION") { 
    ($InputLine,$InputComment)=split(' ',$InputComment,2);

     if ($InputLine eq "UNIFORM") { 
       ampsConfigLib::RedefineMacro("_COMET_GAS_INJECTION_VELOCITY_DIRECTION_MODE_","_COMET_GAS_INJECTION_VELOCITY_DIRECTION_MODE__UNIFORM_","main/Comet.dfn");
     }
  }

  #forces that will be accounted during the simulation
  elsif ($InputLine eq "FORCES") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
    if ($InputLine eq "GRAVITY") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__GRAVITY_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__GRAVITY_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }   
    elsif ($InputLine eq "FRAMEROTATION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__FRAME_ROTATION_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__FRAME_ROTATION_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "RADIATIONPRESSURE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__RADIATION_PRESSURE_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__RADIATION_PRESSURE_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "LORENTZFORCE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__LORENTZ_FORCE_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__LORENTZ_FORCE_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }       
    elsif ($InputLine eq "DRAGFORCE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "DRAGFORCETANGENTIALCOMPONENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE__TANGENTIAL_COMPONENT_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE__TANGENTIAL_COMPONENT_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  
  #dust initial velocity model
  elsif ($InputLine eq "DUSTINITIALVELOCITY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "MODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "CONSTANT") {
        ampsConfigLib::ChangeValueOfVariable("int Comet::DustInitialVelocity::VelocityModelMode","Comet::DustInitialVelocity::Mode::ConstantVelocity","main/Comet.cpp");
      }
      elsif ($InputLine eq "ROTATIONBODY") {
        ampsConfigLib::ChangeValueOfVariable("int Comet::DustInitialVelocity::VelocityModelMode","Comet::DustInitialVelocity::Mode::RotationBody","main/Comet.cpp");        
      }
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "INJECTIONCONSTANTSPEED") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("double Comet::DustInitialVelocity::InjectionConstantVelocity","$InputLine","main/Comet.cpp");   
    } 
    elsif ($InputLine eq "ROTATIONPERIOD") {
      my @all;
      
      $line=$LineOriginal;
      chomp($line);
      $line=~s/[=:!]/ /g;
      
      @all=split(' ',$line);
      ampsConfigLib::ChangeValueOfVariable("double Comet::DustInitialVelocity::RotationPeriod","$all[2]","main/Comet.cpp");    
    } 
    elsif ($InputLine eq "ROTATIONAXIS") {
      my @Axis;
      my @t;
      
      $InputComment=~s/,/ /g;
      @t=split(' ',$InputComment);
      
      #normalize the rotation axis vector
      my $l=sqrt($t[0]*$t[0]+$t[1]*$t[1]+$t[2]*$t[2]);
      
      push(@Axis,$t[0]/$l);
      push(@Axis,$t[1]/$l);
      push(@Axis,$t[2]/$l);
      
      ampsConfigLib::ChangeValueOfArray("double Comet::DustInitialVelocity::RotationAxis\\[3\\]",\@Axis,"main/Comet.cpp");
    }
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }     
  }
  
  #dust charging processes that will be modeled
  elsif ($InputLine eq "DUSTCHARGING") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "ELECTRONCOLLECTIONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ELECTRON_COLLECTION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ELECTRON_COLLECTION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    elsif ($InputLine eq "IONCOLLECTIONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ION_COLLECTION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ION_COLLECTION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    elsif ($InputLine eq "PHOTOELECTRONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    elsif ($InputLine eq "SECONDARYELECTRONEMISSIONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    
    elsif ($InputLine eq "MODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);

      if  ($InputLine eq "OFF")  {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING_MODE_","_DUST__CHARGING_MODE__OFF_","models/dust/Dust.dfn");
        ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__OFF_","pic/picGlobal.dfn");
      }
      elsif  ($InputLine eq "TIMEDEPENDENT") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING_MODE_","_DUST__CHARGING_MODE__TIME_DEPENDENT_","models/dust/Dust.dfn");
        ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "EQUILIBRIUM") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING_MODE_","_DUST__CHARGING_MODE__EQUILIBRIUM_","models/dust/Dust.dfn");
        ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_","pic/picGlobal.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }    
    
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }    
  
  #gravity mode  
  elsif ($InputLine eq "GRAVITY3D") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_3DGRAVITY__MODE_","_3DGRAVITY__MODE__ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_3DGRAVITY__MODE_","_3DGRAVITY__MODE__OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "ACTIVEREGION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_COMET_DUST_USE_ACTIVE_REGION_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_COMET_DUST_USE_ACTIVE_REGION_","_PIC_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "READDRAGFORCEFROMBATL") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_COMET_READ_DRAG_FORCE_FROM_BATL_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_COMET_READ_DRAG_FORCE_FROM_BATL_","_PIC_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "SAVEOUTPUTBINARY") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_SAVING_BINARY_OUTPUT_MODE_","_SAVING_BINARY_OUTPUT_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_SAVING_BINARY_OUTPUT_MODE_","_SAVING_BINARY_OUTPUT_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "READNEUTRALSFROMBINARY") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_READ_NEUTRALS_FROM_BINARY_MODE_","_READ_NEUTRALS_FROM_BINARY_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_READ_NEUTRALS_FROM_BINARY_MODE_","_READ_NEUTRALS_FROM_BINARY_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "NUMBEROFNEUTRALSTOREAD") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("int Comet::CometData::nNeutrals","$InputLine","main/CometData.cpp");
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
  elsif ($InputLine eq "SAMPLEBACKFLUXMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_SAMPLE_BACKFLUX_MODE_","_SAMPLE_BACKFLUX_MODE__ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_SAMPLE_BACKFLUX_MODE_","_SAMPLE_BACKFLUX_MODE__OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "COMPUTEMAXIMUMLIFTABLESIZEMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE_","_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE__ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE_","_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE__OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "NUMERICALLOSSRATEINCREASE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("const double NumericalLossRateIncrease","$InputLine","main/Comet.h");
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
  elsif ($InputLine eq "COMETTEMPERATUREMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "BJORN") {
	  ampsConfigLib::RedefineMacro("_COMET_TEMPERATURE_MODE_","_COMET_TEMPERATURE_MODE__BJORN_","main/Comet.dfn");
      }
      elsif ($InputLine eq "CONSTANT") {
	  ampsConfigLib::RedefineMacro("_COMET_TEMPERATURE_MODE_","_COMET_TEMPERATURE_MODE__CONSTANT_","main/Comet.dfn");
      }
      elsif ($InputLine eq "ANALYTICAL") {
	  ampsConfigLib::RedefineMacro("_COMET_TEMPERATURE_MODE_","_COMET_TEMPERATURE_MODE__ANALYTICAL_ ","main/Comet.dfn");
      }
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
      elsif ($InputLine eq "ANALYTICAL") {
	  ampsConfigLib::RedefineMacro("_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_","_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ANALYTICAL_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "MODELSOURCEDFMS") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_MODEL_SOURCE_DFMS_","_MODEL_SOURCE_DFMS_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_MODEL_SOURCE_DFMS_","_MODEL_SOURCE_DFMS_OFF_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "MULTISPECIESANALYTICALMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_MULTISPECIES_ANALYTICAL_MODE_","_MULTISPECIES_ANALYTICAL_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_MULTISPECIES_ANALYTICAL_MODE_","_MULTISPECIES_ANALYTICAL_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "TRACKINGSURFACEELEMENTMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_TRACKING_SURFACE_ELEMENT_MODE_","_TRACKING_SURFACE_ELEMENT_MODE_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_TRACKING_SURFACE_ELEMENT_MODE_","_TRACKING_SURFACE_ELEMENT_MODE_OFF_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "DFMSRATIOMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_DFMS_RATIO_MODE_","_DFMS_RATIO_MODE_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_DFMS_RATIO_MODE_","_DFMS_RATIO_MODE_OFF_ ","main/Comet.dfn");
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
  elsif ($InputLine eq "DUSTMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__MODE_","_PIC_MODEL__DUST__MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__MODE_","_PIC_MODEL__DUST__MODE__OFF_","pic/picGlobal.dfn");
          ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__OFF_","pic/picGlobal.dfn");
      }
      else {
	  die "Unknown option\n";
      }
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
  
  
  #dust surface injection corection mode
  elsif ($InputLine eq "DUSTINJECTIONCORRECTIONMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
#     print "$InputLine\n";
      
      if (($InputLine eq "INJECTIONRATE") || ($InputLine eq "CORRECTIONFLUXRELATIVEH2O")) {
        if  ($InputLine eq "INJECTIONRATE") {
          ampsConfigLib::RedefineMacro("_DUST_INJECTION_RATE_CORRECTION_MODE_","_DUST_INJECTION_RATE_CORRECTION_MODE__INJECTION_RATE_ ","main/Comet.dfn");
        }
        elsif ($InputLine eq "CORRECTIONFLUXRELATIVEH2O") {
          ampsConfigLib::RedefineMacro("_DUST_INJECTION_RATE_CORRECTION_MODE_","_DUST_INJECTION_RATE_CORRECTION_MODE__RELATIVE_TO_H2O_ ","main/Comet.dfn");
        }
        
        #get the file name 
        my ($s1,$s2,$s3,$s4,$s5);
 
        $LineOriginal=~s/[=():,]/ /g;
        ($s1,$s2,$s3,$s4,$s5)=split(' ',$LineOriginal,5);
        chomp($LineOriginal);
        
        ampsConfigLib::ChangeValueOfVariable("char DustInjectionRateCorrectionFactorFileName\\[\\]","\"".$s4."\"","main/Comet.cpp");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST_INJECTION_RATE_CORRECTION_MODE_","_DUST_INJECTION_RATE_CORRECTION_MODE__OFF_ ","main/Comet.dfn");
      }
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
  }   
  
  #init the location of the Sun from the initial simjulation time string before the beginig of the simulation
  elsif ($InputLine eq "INITSOLARLOCATION") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);

    if  ($InputLine eq "ON") {
      ampsConfigLib::ChangeValueOfVariable("bool Comet::Time::InitSunLocationFlag","true","main/Comet.cpp");
    }
    elsif ($InputLine eq "OFF") {
      ampsConfigLib::ChangeValueOfVariable("bool Comet::Time::InitSunLocationFlag","false","main/Comet.cpp");
    }
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  } 
  
  #recalculate the location of the Sun each iteration during the model run
  elsif ($InputLine eq "RECALCULATESOLARLOCATION") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
     	
    if  ($InputLine eq "ON") {
      ampsConfigLib::ChangeValueOfVariable("bool Comet::Time::RecalculateSunLocationFlag","true","main/Comet.cpp");
    }
    elsif ($InputLine eq "OFF") {
      ampsConfigLib::ChangeValueOfVariable("bool Comet::Time::RecalculateSunLocationFlag","false","main/Comet.cpp");
    }
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  
  #update the boundary condition, based on the  change of solar location.
  elsif ($InputLine eq "UPDATEBOUNDARYCONDITION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON"){
	  ampsConfigLib::ChangeValueOfVariable("bool UpdateBoundaryConditionFlag","true","main/main.cpp");
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
	      
	      if ($s0 eq "UPDATEPERIOD") {
		  ($s0,$line)=split(' ',$line,2);
		  ampsConfigLib::ChangeValueOfVariable("double updateBoundaryConditionPeriod",$s0,"main/main.cpp");
		  
		  $line=~s/(=)/ /;
       	      }
	      else {
		  die "$InputLine: Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
	      }
	  }
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::ChangeValueOfVariable("bool UpdateBoundaryConditionFlag","false","main/main.cpp");
      }
      else {
	  die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
  }
  #starting time of collecting the end of the mission Rosina data 
  elsif ($InputLine eq "SIMULATIONSTARTTIME") {
    
    ($InputLine,$InputComment)=split('!',$line,2);

    $InputLine=~s/[=(),]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);
    
    
    
    ampsConfigLib::ChangeValueOfVariable("char Comet::Time::SimulationStartTimeString\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputComment."\"","main/Comet.cpp");
  }
  
  #starting time of collecting the end of the mission Rosina data 
  elsif ($InputLine eq "STARTSAMPLINGTIME") {
    ($InputLine,$InputComment)=split('!',$line,2);

    $InputLine=~s/[=(),]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);

    ampsConfigLib::ChangeValueOfVariable("char RosinaSample::SamplingTimeInterval::StartSamplingTime\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputComment."\"","main/RosinaMeasurements.cpp");
  }
  
  #end time of collecting the end of the mission Rosina data
  elsif ($InputLine eq "ENDSAMPLINGTIME") {
    ($InputLine,$InputComment)=split('!',$line,2);

    $InputLine=~s/[=(),]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);

    ampsConfigLib::ChangeValueOfVariable("char RosinaSample::SamplingTimeInterval::EndSamplingTime\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputComment."\"","main/RosinaMeasurements.cpp");
  }  
  
  #output the sample length after each data file output 
  elsif ($InputLine eq "INCREASESAMPLINGLENGTH") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "ON" ) {
      ampsConfigLib::ChangeValueOfVariable("bool Comet::IncrementSamplingLengthFlag","true","main/Comet.cpp");     
    }
    elsif ($InputLine eq "OFF" ) {
      ampsConfigLib::ChangeValueOfVariable("bool Comet::IncrementSamplingLengthFlag","false","main/Comet.cpp");     
    }
    else {
     die "Option is unknown (_LINE_), line=$InputFileLineNumber, option=$InputLine ($InputFileName)\n";
    } 
  }
  
  #update the Sun location section
  elsif ($InputLine eq "UPDATESUNLOCATION") {
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
 
    $InputLine=~s/[=]/ /g;
    ($s0,$InputLine)=split(' ',$InputLine,2);
    ($s0,$InputLine)=split(' ',$InputLine,2);
  
    while (defined $InputLine) {     
      ($s0,$InputLine)=split(' ',$InputLine,2);
      
      if ($s0 eq "TIMEINCREMENT") {
        ($s0,$InputLine)=split(' ',$InputLine,2);
        
        ampsConfigLib::ChangeValueOfVariable("double Comet::SunLocationUpdate::TimeIncrement",$s0,"main/SunLocation.cpp");        
      }
      elsif ($s0 eq "FIRSTUPDATEOUTPUTCYCLENUMBER") {
        ($s0,$InputLine)=split(' ',$InputLine,2);
        
        ampsConfigLib::ChangeValueOfVariable("int Comet::SunLocationUpdate::FirstUpdateOutputCycleNumber",$s0,"main/SunLocation.cpp");        
      }
      elsif ($s0 eq "OUTPUTCYCLESTEP") {
        ($s0,$InputLine)=split(' ',$InputLine,2);
        
        ampsConfigLib::ChangeValueOfVariable("int Comet::SunLocationUpdate::OutputCycleStep",$s0,"main/SunLocation.cpp");        
      }
      elsif ($s0 eq "STARTTIME") {
        ($s0,$InputLine)=split(' ',$InputLine,2);
        
        ampsConfigLib::ChangeValueOfVariable("char Comet::SunLocationUpdate::StartTime\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$s0."\"","main/SunLocation.cpp");        
      }
      else {
        die "Option is unknown ( sdfsd __LINE__ ), line=$InputFileLineNumber, option=$InputLine ($InputFileName)\n";
      }          
    }
  }
  
  elsif ($InputLine eq "#ENDBLOCK") {
      last;
  }
   
  else {
    die "Option is unknown, line=$InputFileLineNumber, option=$InputLine ($InputFileName)\n";
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
