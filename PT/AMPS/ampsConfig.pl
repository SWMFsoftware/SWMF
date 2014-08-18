#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#$Id$
 
use strict;
use warnings;

use ampsConfigLib;
use constant {true => 1, false =>0};

my $loadedFlag_MainBlock=0;
my $loadedFlag_SpeciesBlock=0;
my $loadedFlag_BackgroundSpeciesBlock=0;

my $InputFileNameDefault="moon.input"; #"cg.input"; #"mercury.input-test"; #"moon.input"; #"mercury.input"; #"moon.input";
my $InputFileName;


my $line;
my $InputFileLineNumber=0;
my $InputLine;
my $InputComment;
my $s0;
my $s1;
my $s2;
my $FileName;

#Main Block Global Variables
my $TotalSpeciesNumber=0;
my @SpeciesList;

#compile the code
my $CompileProcessedCodeFlag=1;

#Location of the local working vertion of the code that will be compiled 
$ampsConfigLib::WorkingSourceDirectory="src";


#location of the code distribution
my $SourceDirectory=" ";

#the location of the project sprcific sources
my $ProjectSpecificSourceDirectory="main";

#compilation mode: stand along of a part of SWMF
my $CompilationMode="AMPS";

#the number of the compiling threads. The variable is UNDEFINED by default
my $nCompilingThreads;

#markers
my $MARKER__SPECIES_MACRO_DEFINIETION_USED_IN_SIMULATION;#the list of macro used to denote the species used in the simulation (pic.h)

#relate the variables in the input file with the variables and source files of the code
my %InputFileVariable2CodeVariable=("TotalSpeciesNumber"     =>      "static const int nTotalSpecies");

my %InputFileVariable2CodeSourceFile=("TotalSpeciesNumber"     =>    "pic/pic.h");

#read the argument line of the script 
for (my $i=0;$i<$#ARGV + 1;$i++) {
  if ($ARGV[$i] eq "-input") {
    $i++;
    $InputFileName=$ARGV[$i];
  }
  elsif ($ARGV[$i] eq "-mode") {
    $i++;
    $CompilationMode=$ARGV[$i];
  }
  elsif ($ARGV[$i] eq "-no-compile") {
    $CompileProcessedCodeFlag=0;
  }
  elsif ($ARGV[$i] eq "-j") {
    $nCompilingThreads=$ARGV[++$i];
  }
  else {
    die "Unknown argument \"$ARGV[$i]\"\n";
  }
  
}


#if the input file is not determine in the argument line: 1. search for file .amps.InputFileName.txt that containes the inputfile name or 2. use the default input file name
if (!defined($InputFileName)) {
  if (-f ".amps.InputFile.txt") {
    open(FILE,"<.amps.InputFileName.txt") || die "Cannot open .amps.InputFileName.txt\n";
    $InputFileName=<FILE>;
    close(FILE);
  }
  else {
    $InputFileName=$InputFileNameDefault;
  }
}

open(FILE,">.amps.InputFileName.txt");
print FILE "$InputFileName";
close(FILE);
 

#output basic parameters of the code configuration 
print "InputFile: $InputFileName\n"; 

#assemble the input file 
open (AssembledInputFile,">","$InputFileName.Assembled");

AssambleInputFile($InputFileName);
close AssembledInputFile;

open (InputFile,"<","$InputFileName.Assembled") || die "Cannot find file \"$InputFileName.Assembled\"\n";


#read the assembled input file
print "Preprocessing AMPS sources\n";

while ($line=<InputFile>) {
  ($InputFileLineNumber,$FileName)=split(' ',$line);
  $line=<InputFile>;
  
  ($InputLine,$InputComment)=split('!',$line,2);
  $InputLine=uc($InputLine);
  chomp($InputLine);
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;

  #Call subroutines that reads particular blocks of the input file
  if ($InputLine eq "#MAIN") {
    $loadedFlag_MainBlock=1;
    ReadMainBlock();
  }
  elsif ($InputLine eq "#SPECIES") {
    $loadedFlag_SpeciesBlock=1;
    ReadSpeciesBlock();
  }
  elsif ($InputLine eq "#BACKGROUNDSPECIES") {
    $loadedFlag_BackgroundSpeciesBlock=1;
    
    if ($loadedFlag_SpeciesBlock eq 0) {
      die "Species block MUST be loaded before BackgroundSpecies\n";
    }
    
    ReadBackgroundAtmosphereBlock();
  }
  elsif ($InputLine eq "#PARTICLECOLLISIONS") {
    ParticleCollisionModel();
  }
  elsif ($InputLine eq "#IDF") {
    ReadIDF();
  }
  elsif ($InputLine eq "#UNIMOLECULARREACTIONS") {
    ReadUnimolecularReactions();
  }
  elsif ($InputLine eq "#SAMPLING") {
    Sampling();
  }  
  elsif ($InputLine eq "#USERDEFINITIONS") {
    UserDefinitions();
  } 
  elsif ($InputLine eq "#GENERAL") {
    ReadGeneralBlock();
  }  
  elsif ($InputLine eq "#BLOCK") {
    #call a user defined processor of a block in the input file 
    my $BlockProcessor;
      
    chomp($line);
    ($InputLine,$InputComment)=split('!',$line,2);     
    ($InputLine,$BlockProcessor)=split(' ',$InputLine,2);    
    $BlockProcessor=~s/ //g;
    $BlockProcessor=$ampsConfigLib::WorkingSourceDirectory."/".$BlockProcessor;
      
    if (! -e $BlockProcessor) {
      die "Cannot find the input file block processor $BlockProcessor\n";
    }
        
    #extract the block into a new file and call the block processor
    open (BLOCK, ">","$InputFileName.Assembled.Block"); 
      
    while ($line=<InputFile>) {
      $InputFileLineNumber++;
    
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=uc($InputLine);
      chomp($InputLine);    
                 
      print BLOCK "$line";  
      
#      print "$InputFileLineNumber: $InputLine\n";
#      if ($InputFileLineNumber==104) {
#        print "!!!\n";
#      }
        
      if ($InputLine eq "#ENDBLOCK") {
        last;
      }  
    }
      
    close BLOCK;
    
    my $cmd="perl $BlockProcessor $InputFileName.Assembled.Block $ampsConfigLib::WorkingSourceDirectory";
    print "Call: $cmd\n";
  
    system($cmd) and do {
      die "Execution of $BlockProcessor was not succesful\n";
    };
  }
  
  elsif ($InputLine eq "#END") {
    last;
  }
  else {
    die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
  }
  
}

#read settings from .ampsConfig
if (-e ".ampsConfig.Settings") {
  ampsConfigSettings();
}

#modify the makefile for the compilation mode
if ($CompileProcessedCodeFlag==1) {
  my @makefilelines;
        
  open (MAKEFILE,"<Makefile.local") || die "Cannot open Makefile.local\n";
  @makefilelines=<MAKEFILE>;
  close (MAKEFILEFILE);
        
  foreach (@makefilelines) {
      if ($_=~/\${CC} -o \${EXE}/) {
        $_=~s/^#/\t/;
        
        if ($CompilationMode ne "AMPS") {
          $_=~s/^\t/#/;
        }
      }
  }
        
  open (MAKEFILEFILE,">Makefile.local") || die "Cannot open Makefile.local\n";
  print MAKEFILEFILE @makefilelines;
  close (MAKEFILEFILE);

  #compile the code 
  print "Compile the code\n";
  
  if (defined $nCompilingThreads) {
    system("make -j $nCompilingThreads"); 
  }
  else {
    system("make");
  }
}


#=============================== INCLUDE an additional input file to the assemble =============================
sub AssambleInputFile {
  my $line;
  my $InputLine;
  my $s0;
  
  my $LineBreak=0,
  my $nline=0;
  
  my $StatmentEnd=1;
  
  open (my $Input,"<",$_[0]) || die "Cannot find file \"$_[0]\"\n";
  
  while ($line=<$Input>) {
   $nline++;
        
   ($line,$InputLine)=split('!',$line,2); 
   $line=~s/^\s+//; #remove spaces from the begining of the line
   $line=~s/\s+$//; #remove spaces from the end of the line
   
   if ((!defined $line) || ($line eq '')) {
     if ($LineBreak==1) {
       print "The next line after a line break cannot be empty or fully commented (file=$_[0],line=$nline)\n";
       die;
     }
   }
   
   $LineBreak=0;
   if ($line=~/\\/) {
     $LineBreak=1;
   }
    
   $InputLine=uc($line);
   chomp($InputLine);
   ($s0,$InputLine)=split(' ',$InputLine,2);
   
   if (!defined $s0) {
     $s0="nothing";
   }
   
   if ($s0 eq "#INCLUDE") {
     #include a file into the main input file
     chomp($line);
     $line=~s/\s+$//;
     ($s0,$line)=split(' ',$line,2);
     ($s0,$line)=split(' ',$line,2);     
     
     AssambleInputFile($s0);
   }
   else {   
     no warnings 'uninitialized';
       
     if (length $line) {
       if ($StatmentEnd==1) {
         print AssembledInputFile "$nline $_[0]\n";
       }
       
       if ($line=~m/\\\\$/) {
         $line=~s/\\\\$//;
         chomp($line);
         print AssembledInputFile "$line  ";
         $StatmentEnd=0;
       }
       else {
         print AssembledInputFile "$line\n";
         $StatmentEnd=1;
       }
     }
   }
}
  
}


#=============================== Read Main Block =============================
sub ReadMainBlock {
  my $s0;
  my $s1;
    
  #set up the prefix,error log,diagnistic streem
  my $Prefix='';
  my $ErrorLog='error.log';
  my $DiagnosticStream='stdout';
  my $OutputDirectory='.';
  my $StdoutErrorLog='_STDOUT_ERRORLOG_MODE__ON_';
  
  my $SimulationTimeStepMode; #='_SINGLE_GLOBAL_TIME_STEP_';
  my $SimulationParticleWeightMode; #='_SINGLE_GLOBAL_PARTICLE_WEIGHT_';
  my $SimulationParticleWeightCorrectionMode; #='_INDIVIDUAL_PARTICLE_WEIGHT_OFF_';
  
  my $CouplingMode;
  my $TrajectoryIntegrationCheckBlockFaceIntersection;
  
  #force the repeatable execution path
  my $ForceRepatableExecutionPath=0;
  
  #debbuger mode of the execution 
  my $DebuggerMode=0;
    
  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[=,]/ /g;
    
    #get the first word in the sequence
    ($s0,$s1)=split(' ',$InputLine,2);
    
    if ($s0 eq "SPECIESLIST") {
      $TotalSpeciesNumber=0;
      @SpeciesList=split(' ',$s1);
      $MARKER__SPECIES_MACRO_DEFINIETION_USED_IN_SIMULATION=" ";
      
      foreach (@SpeciesList) {
        $MARKER__SPECIES_MACRO_DEFINIETION_USED_IN_SIMULATION=$MARKER__SPECIES_MACRO_DEFINIETION_USED_IN_SIMULATION."\n#undef _".$_."_SPEC_\n#define _".$_."_SPEC_ $TotalSpeciesNumber\n";
        $TotalSpeciesNumber++;
      }        
    }
    elsif ($s0 eq "MAKEFILE") {
      #substitute the following line in the Makefile.local 
      my $Variable;
      my $Argument;
      my @MakeFileContent;
      
      my $FoundVariable=0;
      
      ($line,$Argument)=split('!',$line,2);
      $Variable=$line;
      $Variable=~s/=/ /;
      ($Variable,$Argument)=split(' ',$Variable,2);
      ($Variable,$Argument)=split(' ',$Argument,2);
      
      if (-e "Makefile.local") {     
        open (MAKEFILE,"<","Makefile.local") || die "Cannot open Makefile.local\n";
        @MakeFileContent=<MAKEFILE>;
        close(MAKEFILE);
        
        my $l0=" ";
        my $l1=" ";
        
        foreach (@MakeFileContent) {
          $l1=$_;
          $l1=~s/ //g;
          $l1=~s/=/ /;
          
          ($l0,$l1)=split(' ',$l1,2);
          
          if (!defined $l0) {
            $l0="nothing";
          }
          
          if ($l0 eq $Variable) {
            ($l0,$l1)=split(' ',$line,2);
            $_=$l1;
            $FoundVariable=1;
          }
          #elsif ($l0 eq "CWD") {
          #  $l1=`pwd`;
          #  $_="CWD=$l1";
          #}
          
        }
        
        if ($FoundVariable == 0) {
          #the definition of the variable is not found -> add it to Makefile.local
          push(@MakeFileContent,$Variable."=".$Argument);
        }
      }
      
      #write the updated copy of the Makefile.local 
      open (MAKEFILE,">","Makefile.local");

      foreach (@MakeFileContent) {
        print MAKEFILE "$_";
      }
      
      close (MAKEFILE);
    }
    
    elsif ($s0 eq "FORCEREPEATABLESIMULATIONPATH") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        $ForceRepatableExecutionPath=1;
      }
      elsif ($s0 eq "OFF") {
        $ForceRepatableExecutionPath=0;
      }
      else {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    elsif ($s0 eq "DEBUGGERMODE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        $DebuggerMode=1;
      }
      elsif ($s0 eq "OFF") {
        $DebuggerMode=0;
      }
      else {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    
    
    
        
    elsif ($s0 eq "COUPLERMODE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "OFF") {$CouplingMode="_PIC_COUPLER_MODE__OFF_";}
      elsif ($s0 eq "ICES") {$CouplingMode="_PIC_COUPLER_MODE__ICES_";}
      elsif ($s0 eq "SWMF") {$CouplingMode="_PIC_COUPLER_MODE__SWMF_";}
      else {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    
    elsif ($s0 eq "WORKINGSOURCEDIRECTORY") {
      chomp($line);
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=~s/ //g;
      $InputLine=~s/=/ /;
      
      ($InputLine,$ampsConfigLib::WorkingSourceDirectory)=split(' ',$InputLine,2);
    }
    elsif ($s0 eq "SOURCEDIRECTORY") {
      chomp($line);
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=~s/ //g;
      $InputLine=~s/=/ /;
      
      ($InputLine,$SourceDirectory)=split(' ',$InputLine,2);
    } 
    
    
    elsif ($s0 eq "PROJECTSOURCEDIRECTORY") {
      chomp($line);
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=~s/ //g;
      $InputLine=~s/=/ /;
      
      ($InputLine,$ProjectSpecificSourceDirectory)=split(' ',$InputLine,2);
    }
       
    
    elsif ($s0 eq "TRAJECTORYINTERSECTIONWITHBLOCKFACES") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if(    $s0 eq "ON") {
	  $TrajectoryIntegrationCheckBlockFaceIntersection="_PIC_MODE_ON_";
      }elsif($s0 eq "OFF") {
	  $TrajectoryIntegrationCheckBlockFaceIntersection="_PIC_MODE_OFF_";
      }else{
	  die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }       
    
    elsif ($s0 eq "STDOUTERRORLOG") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        $StdoutErrorLog='_STDOUT_ERRORLOG_MODE__ON_';
      }
      elsif ($s0 eq "OFF") {
        $StdoutErrorLog='_STDOUT_ERRORLOG_MODE__OFF_';
      }
      else {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }    
    
    
    elsif ($s0 eq "TIMESTEPMODE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "SINGLEGLOBALTIMESTEP") {         
        $SimulationTimeStepMode='_SINGLE_GLOBAL_TIME_STEP_';
      }
      elsif ($s0 eq "SPECIESGLOBALTIMESTEP") {
        $SimulationTimeStepMode='_SPECIES_DEPENDENT_GLOBAL_TIME_STEP_';
      }
      elsif ($s0 eq "SINGLELOCALTIMESTEP") {
        $SimulationTimeStepMode='_SINGLE_LOCAL_TIME_STEP_';
      }
      elsif ($s0 eq "SPECIESLOCALTIMESTEP") {
        $SimulationTimeStepMode='_SPECIES_DEPENDENT_LOCAL_TIME_STEP_';
      }
      else {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }    
    elsif ($s0 eq "PARTICLEWEIGHTMODE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "SINGLEGLOBALPARTICLEWEIGHT") {
        $SimulationParticleWeightMode='_SINGLE_GLOBAL_PARTICLE_WEIGHT_';
      }
      elsif ($s0 eq "SPECIESGLOBALPARTICLEWEIGHT") {
        $SimulationParticleWeightMode='_SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_';
      }
      elsif ($s0 eq "SINGLELOCALPARTICLEWEIGHT") {
        $SimulationParticleWeightMode='_SINGLE_LOCAL_PARTICLE_WEIGHT_';
      }
      elsif ($s0 eq "SPECIESLOCALPARTICLEWEIGHT") {
        $SimulationParticleWeightMode='_SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_';
      }
      else {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    elsif ($s0 eq "PARTICLEWEIGHTCORRECTIONMODE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        $SimulationParticleWeightCorrectionMode='_INDIVIDUAL_PARTICLE_WEIGHT_ON_'
      }
      elsif ($s0 eq "OFF") {
        $SimulationParticleWeightCorrectionMode='_INDIVIDUAL_PARTICLE_WEIGHT_OFF_'
      }
      else {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
   
    elsif ($s0 eq "ERRORLOG") {
      chomp($line);
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=~s/ //g;
      $InputLine=~s/=/ /;
      
      ($InputLine,$ErrorLog)=split(' ',$InputLine,2);
    } 
    elsif ($s0 eq "PREFIX") {
      chomp($line);
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=~s/ //g;
      $InputLine=~s/=/ /;
      
      ($InputLine,$Prefix)=split(' ',$InputLine,2);
      $Prefix=$Prefix.":";
    }
    elsif ($s0 eq "DIAGNOSTICSTREAM") {
      chomp($line);
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=~s/ //g;
      $InputLine=~s/=/ /;
      
      ($InputLine,$DiagnosticStream)=split(' ',$InputLine,2);
      
      if ($DiagnosticStream eq 'screen') {
        $DiagnosticStream='stdout';
      }
    }  
    elsif ($s0 eq "OUTPUTDIRECTORY") {
      chomp($line);
      ($InputLine,$InputComment)=split('!',$line,2);
      $InputLine=~s/ //g;
      $InputLine=~s/=/ /;
      
      ($InputLine,$OutputDirectory)=split(' ',$InputLine,2);
    }  
    
    elsif ($s0 eq "#ENDMAIN") {
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
  }
  
  
  #set working directory in the makefile 
  if (-e "Makefile.local") {
    my @MakeFileContent;
    my $FoundWSD=0;
    
    open (MAKEFILE,"<","Makefile.local") || die "Cannot open Makefile.local\n";
    @MakeFileContent=<MAKEFILE>;
    close(MAKEFILE);
    
    open (MAKEFILE,">","Makefile.local");

    foreach (@MakeFileContent) {
      if ($_=~m/WSD=/) {
        $_="WSD=".$ampsConfigLib::WorkingSourceDirectory."\n";
        $FoundWSD=1;
      }
      
      print MAKEFILE "$_";
    }
      
    close (MAKEFILE);
    
    #check if the WSD is found
    if ($FoundWSD eq 0 ) {
      open (MAKEFILE,"<","Makefile.local") || die "Cannot open Makefile.local\n";
      @MakeFileContent=<MAKEFILE>;
      close(MAKEFILE);
    
      open (MAKEFILE,">","Makefile.local");   
      print MAKEFILE "WSD=$ampsConfigLib::WorkingSourceDirectory\n";  
    
      foreach (@MakeFileContent) {
        print MAKEFILE "$_";      
      }  
      
      
      close (MAKEFILE);
    }
  }
  
  
  #remove the temporary source directory and copy there a fresh copy of the code distribution
  if ( -d $ampsConfigLib::WorkingSourceDirectory ) {
    `rm -f -r $ampsConfigLib::WorkingSourceDirectory`;
  }
  
  `cp -r $SourceDirectory $ampsConfigLib::WorkingSourceDirectory`;
 
  if ( -d $ProjectSpecificSourceDirectory ) {
    `cp -r $ProjectSpecificSourceDirectory $ampsConfigLib::WorkingSourceDirectory/main`;
  }
  
  #setup the time-step and particle-weight modes
  if (defined $SimulationTimeStepMode) {ampsConfigLib::RedefineMacro("_SIMULATION_TIME_STEP_MODE_",$SimulationTimeStepMode,"pic/picGlobal.dfn");}
  if (defined $SimulationParticleWeightMode) {ampsConfigLib::RedefineMacro("_SIMULATION_PARTICLE_WEIGHT_MODE_",$SimulationParticleWeightMode,"pic/picGlobal.dfn");}
  if (defined $SimulationParticleWeightCorrectionMode) {ampsConfigLib::RedefineMacro("_INDIVIDUAL_PARTICLE_WEIGHT_MODE_",$SimulationParticleWeightCorrectionMode,"pic/picGlobal.dfn");}
  
  #check intersection of particle trajectory with the boundaries of the blocks (needed with local timestepping, local weight and when the internal surface has fine fiatures)
  if (defined $TrajectoryIntegrationCheckBlockFaceIntersection) {ampsConfigLib::RedefineMacro("_PIC__PARTICLE_MOVER__CHECK_BLOCK_FACE_INTERSECTION__MODE_",$TrajectoryIntegrationCheckBlockFaceIntersection,"pic/picGlobal.dfn");}
  
  
  #change prefix,error log file and diagnostic stream
  ampsConfigLib::RecursiveSubstitute('\$PREFIX:',$Prefix,$ampsConfigLib::WorkingSourceDirectory);
  ampsConfigLib::RecursiveSubstitute('\$ERRORLOG',$ErrorLog,$ampsConfigLib::WorkingSourceDirectory);
  ampsConfigLib::ChangeValueOfVariable("char PIC::DiagnospticMessageStreamName\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$DiagnosticStream."\"","pic/pic_init_const.cpp");
  ampsConfigLib::ChangeValueOfVariable("char PIC::OutputDataFileDirectory\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$OutputDirectory."\"","pic/pic_init_const.cpp");
  
  #redefine the value for the macro controlling output error messages into stdout
  ampsConfigLib::RedefineMacro("_STDOUT_ERRORLOG_MODE_",$StdoutErrorLog,"general/specfunc.h");
  
  #update variables of the code distribution with the parameters of the "main block" of the input file 
  ampsConfigLib::ChangeValueOfVariable($InputFileVariable2CodeVariable{"TotalSpeciesNumber"},$TotalSpeciesNumber,$InputFileVariable2CodeSourceFile{"TotalSpeciesNumber"});
  
  #redefine the value of the macro that allows loading the user-defined table of the macros that describe the species used in the simulation 
  ampsConfigLib::RedefineMacro("_PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE_","_PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE__ON_","pic/picGlobal.dfn");
  
  #redefine the value of the macro that determined the coupling of AMPS
  if (defined $CouplingMode) {ampsConfigLib::RedefineMacro("_PIC_COUPLER_MODE_",$CouplingMode,"pic/picGlobal.dfn");}


  #add markers to the code
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/pic/pic.h") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/pic/pic.h\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/pic/pic.h") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/pic/pic.h\n";
  
  foreach (@FileContent) {
    $_=~s/\$MARKER:SPECIES-MACRO-DEFINIETION-USED-IN-SIMULATION\$/$MARKER__SPECIES_MACRO_DEFINIETION_USED_IN_SIMULATION/;
      
    print FILEOUT "$_";
  }
   
  close (FILEOUT);
  
  
  #set the debugger mode
  if ($DebuggerMode == 0) {
    ampsConfigLib::RedefineMacro("_PIC_DEBUGGER_MODE_","_PIC_DEBUGGER_MODE_OFF_","pic/picGlobal.dfn");
  }
  else {
    ampsConfigLib::RedefineMacro("_PIC_DEBUGGER_MODE_","_PIC_DEBUGGER_MODE_ON_","pic/picGlobal.dfn");
  }
  
  #set the conditions for a repeables execution path of the code
  if ($ForceRepatableExecutionPath == 1) {
    ampsConfigLib::RedefineMacro("_PIC_DEBUGGER_MODE_","_PIC_DEBUGGER_MODE_ON_","pic/picGlobal.dfn");
    ampsConfigLib::RedefineMacro("_PIC_DYNAMIC_LOAD_BALANCING_MODE_","_PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_","pic/picGlobal.dfn");
  }

}

#=============================== User Definitions ==================
sub UserDefinitions {

  print "ampsConfig.pl: Warning: the user definition section is obsolete and has no effect on neither configuration nor execution of AMPS\n";

  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[=,]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);
    
    if ($InputLine eq "PIC") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC_USER_DEFINITION_MODE_","_PIC_USER_DEFINITION_MODE__ENABLED_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_USER_DEFINITION_MODE_","_PIC_USER_DEFINITION_MODE__DISABLED_","pic/picGlobal.dfn");
      }
      else {
        die "The option is unknown\n";
      }
    }
    elsif ($InputLine eq "MESH") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_AMR__LOAD_USER_DEFINITION__MODE_","_AMR__LOAD_USER_DEFINITION__MODE__ON_","meshAMR/meshAMRdef.h");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_AMR__LOAD_USER_DEFINITION__MODE_","_AMR__LOAD_USER_DEFINITION__MODE__OFF_","meshAMR/meshAMRdef.h");
      }
      else {
        die "The option is unknown\n";
      }
    }
      
      
    elsif ($InputLine eq "#ENDUSERDEFINITIONS") {
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    
  }
}


#=============================== General Block ==================
sub ReadGeneralBlock {

  
  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[=,]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);
    
   
       
    
    if ($InputLine eq "MAXMESHREFINMENTLEVEL") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::RedefineMacro("_MAX_REFINMENT_LEVEL_",$InputLine,"meshAMR/meshAMRdef.h");
    }
    elsif ($InputLine eq "REFERENCEINJECTIONPARTICLENUMBER") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("double PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber",$InputLine,"pic/pic_weight_time.cpp");
    }
    
    elsif ($InputLine eq "CUTCELLVOLUMECALCULATIONMAXREFINMENTLEVEL") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::RedefineMacro("_AMR__CUT_CELL_VOLUME_CALCULATION__MAX_REFINMENT_LEVEL_",$InputLine,"meshAMR/meshAMRdef.h");
    }

    elsif ($InputLine eq "TESTRUNTOTALITERATION") {
	($InputLine,$InputComment)=split(' ',$InputComment,2);
	ampsConfigLib::ChangeValueOfVariable("int PIC::ModelTestRun::nTotalIteraction",$InputLine,"pic/pic_init_const.cpp");
    }

    elsif ($InputLine eq "INITIALSAMPLELENGTH") {
	($InputLine,$InputComment)=split('!',$line,2);
	chomp($InputLine);
	$InputLine=~s/[=();]/ /g;

	($s0,$s1,$s2)=split(' ',$InputLine,3);
	ampsConfigLib::ChangeValueOfVariable("long int PIC::RequiredSampleLength",$s1,"pic/pic.cpp");
    }

    elsif ($InputLine eq "BLOCKCELLS") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/[,=();]/ /g;

      ($InputLine,$s0,$s1,$s2,$InputComment)=split(' ',$InputLine,5);
      ampsConfigLib::RedefineMacro("_BLOCK_CELLS_X_",$s0,"meshAMR/meshAMRdef.h");
      ampsConfigLib::RedefineMacro("_BLOCK_CELLS_Y_",$s1,"meshAMR/meshAMRdef.h");
      ampsConfigLib::RedefineMacro("_BLOCK_CELLS_Z_",$s2,"meshAMR/meshAMRdef.h");
    } 
    elsif ($InputLine eq "GHOSTCELLS") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/[,=();]/ /g;

      ($InputLine,$s0,$s1,$s2,$InputComment)=split(' ',$InputLine,5);
      ampsConfigLib::RedefineMacro("_GHOST_CELLS_X_",$s0,"meshAMR/meshAMRdef.h");
      ampsConfigLib::RedefineMacro("_GHOST_CELLS_Y_",$s1,"meshAMR/meshAMRdef.h");
      ampsConfigLib::RedefineMacro("_GHOST_CELLS_Z_",$s2,"meshAMR/meshAMRdef.h");
    }

    elsif ($InputLine eq "ENFORCEREQUESTEDMESHRESOLUTION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_AMR_ENFORCE_CELL_RESOLUTION_MODE_","_AMR_ENFORCE_CELL_RESOLUTION_MODE_ON_","meshAMR/meshAMRdef.h");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_AMR_ENFORCE_CELL_RESOLUTION_MODE_","_AMR_ENFORCE_CELL_RESOLUTION_MODE_OFF_","meshAMR/meshAMRdef.h");
      }
      else {
        die "The option is unknown\n";
      }
    }
  
    elsif ($InputLine eq "CONTROLPARTICLEINSIDENASTRANSURFACE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC_CONTROL_PARTICLE_INSIDE_NASTRAN_SURFACE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_CONTROL_PARTICLE_INSIDE_NASTRAN_SURFACE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }
      else {
        die "The option is unknown\n";
      }
    }  
  
    elsif ($InputLine eq "NASTRANSURFACEUSERDATA") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_","_USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_ON_","meshAMR/meshAMRdef.h");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_","_USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_OFF_","meshAMR/meshAMRdef.h");
      }
      else {
        die "The option is unknown\n";
      }
    }      
      
    elsif ($InputLine eq "#ENDGENERAL") {
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    
  }
}

#=============================== Sample particle data ==================
sub Sampling {

  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[=,]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);
  
  
    if ($InputLine eq "SAMPLEPARALLELTANGENTIALKINETICTEMPERATURE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
            
      
      if ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_","_PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "ON") {
        my $t;
        
        $InputComment=~s/[()]/ /g;
        ($t,$InputLine,$InputComment)=split(' ',$InputComment,3);
        
        if ($InputLine eq "CONST") {
          my @x;
          
          @x=split(' ',$InputComment);
          ampsConfigLib::RedefineMacro("_PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_","_PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__CONSTANT_DIRECTION_ORIGIN_","pic/picGlobal.dfn");
          ampsConfigLib::ChangeValueOfArray("static const double constNormalDirection__SampleParallelTangentialTemperature\\[3\\]",\@x,"pic/pic.h");
          
        }
        elsif ($InputLine eq "FUNCTION") {
          die "function is not implemented\n";
        }
        else {
          die "Cannot recognize the option\n";
        }
        
        
#       ampsConfigLib::RedefineMacro("_PIC__PARTICLE_COLLISION_MODEL__MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }  
      else {
        die "Unknown option\n";
      }
    }
      
    elsif ($InputLine eq "#ENDSAMPLING") {
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    
  }
}


#=============================== Read Particle Collision Model ==================
sub ParticleCollisionModel {
  my $ModelIsOnFlag=1;

  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[=,]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);
  
  
    if ($InputLine eq "MODEL") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "HS") {
        ampsConfigLib::RedefineMacro("_PIC__PARTICLE_COLLISION_MODEL__MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
        ampsConfigLib::RedefineMacro("_PIC__PARTICLE_COLLISION_MODEL_","_PIC__PARTICLE_COLLISION_MODEL__HS_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC__PARTICLE_COLLISION_MODEL__MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
        $ModelIsOnFlag=0;
      }  
      else {
        die "Unknown option\n";
      }
    }
    
    elsif ($InputLine eq "SAMPLECOLLISIONFREQUENCY") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__","_PIC_MODE_ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }
      else {
        die "Cannot recognize the option\n";
      }
    }
  
    elsif ($InputLine eq "COLLISIONCROSSSECTION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "FUNCTION") {
        die "not implemented\n";
      }
      elsif ($InputLine eq "CONST") {
        my @CrossSectionTable;
        my ($t0,$t1,$t2,$t3,$s0,$s1,$cnt);
        
        $InputComment=~s/[(),=]/ /g;
                
        for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
          my @tmp=((0)x$TotalSpeciesNumber);
          push(@CrossSectionTable, [@tmp]);
        }      
  
        while (defined $InputComment) {
          ($t0,$t1,$t2,$t3,$InputComment)=split(' ',$InputComment,5);
          
          $s0=-1;
          $s1=-1;
          $cnt=0;
          
  
          foreach (@SpeciesList) {
            if ($_ eq $t1) {
              $s0=$cnt;
            }
            
            if ($_ eq $t2) {
              $s1=$cnt;
            }
                
            $cnt++;
          }
          
          if ($ModelIsOnFlag == 1) {    
            if (($s0 eq -1) || ($s1 eq -1)) {
              die "Cannot recognize model specie (line=$InputFileLineNumber)\n";
            }
                
            $CrossSectionTable[$s0][$s1]=$t3;
            $CrossSectionTable[$s1][$s0]=$t3;    
          }              
        }
        
        #create the new array defienition and substitute the default value
        my $newCrossSectionTable="const double ConstantCollisionCrossSectionTable[PIC::nTotalSpecies][PIC::nTotalSpecies]={";
        
        for ($s0=0;$s0<$TotalSpeciesNumber;$s0++) {
          $newCrossSectionTable=$newCrossSectionTable."{";
          
          for ($s1=0;$s1<$TotalSpeciesNumber;$s1++) {
            $newCrossSectionTable=$newCrossSectionTable.$CrossSectionTable[$s0][$s1];
            
            if ($s1 != $TotalSpeciesNumber-1) {
              $newCrossSectionTable=$newCrossSectionTable.",";
            }
            else {
              $newCrossSectionTable=$newCrossSectionTable."}"
            }
          }
          
          if ($s0 != $TotalSpeciesNumber-1) {
            $newCrossSectionTable=$newCrossSectionTable.",";
          }
          else {
            $newCrossSectionTable=$newCrossSectionTable."};"
          }       
        }
          
          
        ampsConfigLib::SubstituteCodeLine("const double ConstantCollisionCrossSectionTable",$newCrossSectionTable,"pic/pic.h");  
          
        
      }
      else {
        die "Unknown option\n";
      }
    }
    
    elsif ($InputLine eq "#ENDPARTICLECOLLISIONS") {
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    
  }
}

#=============================== Read Internal Degrees of Freedom ==================
sub ReadIDF {
  my $ModelIsOnFlag=1;
  
  my @nTotalVibModes=(0)x$TotalSpeciesNumber;
  my @nTotalRotModes=(0)x$TotalSpeciesNumber;
  my @CharacteristicVibTemp=(0)x$TotalSpeciesNumber;
  my @RotationalZnumber=(0)x$TotalSpeciesNumber;
  
  #temeprature index
  my @TemperatureIndexTable;
  my ($s0,$s1);
  
  for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
    my @tmp=((0)x$TotalSpeciesNumber);
    push(@TemperatureIndexTable, [@tmp]);
  }

  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[()=,]/ /g;
    ($InputLine,$InputComment)=split(' ',$InputLine,2);
  
  
    if ($InputLine eq "MODEL") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "LB") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_","_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "QLB") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_","_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__QLB_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
        $ModelIsOnFlag=0;
      }  
      else {
        die "Unknown option\n";
      }
    }
    
    elsif ($InputLine eq "VVRELAXATION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM__VV_RELAXATION_MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM__VV_RELAXATION_MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }
      else {
        die "Cannot recognize the option\n";
      }
    }
    elsif ($InputLine eq "VTRELAXATION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM__VT_RELAXATION_MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM__VT_RELAXATION_MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }
      else {
        die "Cannot recognize the option\n";
      }
    }
    elsif ($InputLine eq "RTRELAXATION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM__RT_RELAXATION_MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_INTERNAL_DEGREES_OF_FREEDOM__RT_RELAXATION_MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }
      else {
        die "Cannot recognize the option\n";
      }
    }
    
    elsif ($InputLine eq "NVIBMODES") {
      my ($s0,$s1,$nspec);
      
      if ($ModelIsOnFlag == 1) {
        while (defined $InputComment) {      
          ($s0,$s1,$InputComment)=split(' ',$InputComment,3);
             
          if ((defined $s0)&&(defined $s1)) {
            $nspec=ampsConfigLib::GetElementNumber($s0,\@SpeciesList);
            if ($nspec==-1) {
              die "Cannot recognize species '$s0' in line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
            }
          
            $nTotalVibModes[$nspec]=$s1;      
          }
        }
      }
    }
    elsif ($InputLine eq "NROTMODES") {
      my ($s0,$s1,$nspec);
      
      if ($ModelIsOnFlag == 1) {
        while (defined $InputComment) {      
          ($s0,$s1,$InputComment)=split(' ',$InputComment,3);
        
          if ((defined $s0)&&(defined $s1)) {
            $nspec=ampsConfigLib::GetElementNumber($s0,\@SpeciesList);
            if ($nspec==-1) {
              die "Cannot recognize species '$s0' in line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
            }
          
            $nTotalRotModes[$nspec]=$s1;     
          } 
        }     
      } 
    }
    elsif ($InputLine eq "VIBTEMP") {
      my ($s0,$s1,$nspec);
      
      if ($ModelIsOnFlag == 1) {
        while (defined $InputComment) {      
          ($s0,$s1,$InputComment)=split(' ',$InputComment,3);
        
          if ((defined $s0)&&(defined $s1)) {
            $nspec=ampsConfigLib::GetElementNumber($s0,\@SpeciesList);
            if ($nspec==-1) {
              die "Cannot recognize species '$s0' in line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
            }
          
            $CharacteristicVibTemp[$nspec]=$s1;     
          } 
        }      
      }     
    }
    elsif ($InputLine eq "ROTZNUM") {
      my ($s0,$s1,$nspec);
      
      if ($ModelIsOnFlag == 1) {
        while (defined $InputComment) {      
          ($s0,$s1,$InputComment)=split(' ',$InputComment,3);
        
          if ((defined $s0)&&(defined $s1)) {
            $nspec=ampsConfigLib::GetElementNumber($s0,\@SpeciesList);
            if ($nspec==-1) {
              die "Cannot recognize species '$s0' in line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
            }
          
            $RotationalZnumber[$nspec]=$s1;  
          }    
        }   
      }        
    }    
  
    elsif ($InputLine eq "TEMPERATUREINDEX") {     
      if ($ModelIsOnFlag == 1) {     
        ($InputLine,$InputComment)=split(' ',$InputComment,2);     
        $s0=ampsConfigLib::GetElementNumber($InputLine,\@SpeciesList);
        if ($s0==-1) {
          die "Cannot recognize species '$s0' in line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
        }
        
        ($InputLine,$InputComment)=split(' ',$InputComment,2);     
        $s1=ampsConfigLib::GetElementNumber($InputLine,\@SpeciesList);
        if ($s1==-1) {
          die "Cannot recognize species '$s0' in line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
        }      
        
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        $TemperatureIndexTable[$s0][$s1]=$InputLine;
        $TemperatureIndexTable[$s1][$s0]=$InputLine;  
      }
    }

    elsif ($InputLine eq "#ENDIDF") {
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
    
  }
  
  #update the model parameters in the code
  my $nSpeciesMaxVibrationalModes=0;
  
  foreach (@nTotalVibModes) {
    if ($_ > $nSpeciesMaxVibrationalModes) {
      $nSpeciesMaxVibrationalModes=$_;
    }
  }
  
  #create the new TemperatureIndex Table
  my $newTemperatureIndexTable="static const double TemepratureIndex[PIC::nTotalSpecies][PIC::nTotalSpecies]={";
         
  for ($s0=0;$s0<$TotalSpeciesNumber;$s0++) {
    $newTemperatureIndexTable=$newTemperatureIndexTable."{";

    for ($s1=0;$s1<$TotalSpeciesNumber;$s1++) {
      $newTemperatureIndexTable=$newTemperatureIndexTable.$TemperatureIndexTable[$s0][$s1];

      if ($s1 != $TotalSpeciesNumber-1) {
        $newTemperatureIndexTable=$newTemperatureIndexTable.",";
      }
      else {
        $newTemperatureIndexTable=$newTemperatureIndexTable."}"
      }
    }

    if ($s0 != $TotalSpeciesNumber-1) {
      $newTemperatureIndexTable=$newTemperatureIndexTable.",";
    }
    else {
      $newTemperatureIndexTable=$newTemperatureIndexTable."};"
    }     
  }
  
  ampsConfigLib::SubstituteCodeLine("static const double TemepratureIndex",$newTemperatureIndexTable,"pic/pic.h");
  
  
  ampsConfigLib::ChangeValueOfArray("static const int nTotalVibtationalModes\\[\\]",\@nTotalVibModes,"pic/pic.h");
  ampsConfigLib::ChangeValueOfArray("static const int nTotalRotationalModes\\[\\]",\@nTotalRotModes,"pic/pic.h");
  ampsConfigLib::ChangeValueOfArray("static const double CharacteristicVibrationalTemperature\\[\\]",\@CharacteristicVibTemp,"pic/pic.h");
  ampsConfigLib::ChangeValueOfArray("static const double RotationZnumber\\[\\]",\@RotationalZnumber,"pic/pic.h");
  
  ampsConfigLib::ChangeValueOfVariable("static const int nSpeciesMaxVibrationalModes",$nSpeciesMaxVibrationalModes,"pic/pic.h");
  
}

#=============================== Read Background Atmosphere Block ===============
sub ReadBackgroundAtmosphereBlock {
  my $nBackgroundSpecies=0;
  my @BackgroundSpeciesList;
  
  my $s0;
  my $s1;
    
  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[=,]/ /g;
    
    #get the first word in the sequence
    ($s0,$s1)=split(' ',$InputLine,2);
    
    
    if ($s0 eq "BACKGROUNDSPECIES") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_","pic/picGlobal.dfn");
      }
      else {
        die "Cannot recognize the option (line=$InputLine, nline=$InputFileLineNumber)\n";
      }
    }
    
    elsif ($s0 eq "SPECIESLIST") {      
      while (defined $s1) {
        ($s0,$s1)=split(' ',$s1,2);
        push(@BackgroundSpeciesList,$s0);
        $nBackgroundSpecies++;
      }
      
      #create the list of the background species
      open (BACKGROUNDSPECIES,">$InputFileName.Assembled.BackgroundSpecies") || die "Cannot open file $InputFileName.Assembled.BackgroundSpecies\n";
      print BACKGROUNDSPECIES "$nBackgroundSpecies\n";
      
      foreach (@BackgroundSpeciesList) {
        print BACKGROUNDSPECIES "$_\n";
      }
      
      close (BACKGROUNDSPECIES);
      
      #add definitions of the background species to the code
      my $BackgroundSpeciesDefinitionString;
      my $cnt=0;
            
      foreach (@BackgroundSpeciesList) {
        $BackgroundSpeciesDefinitionString=$BackgroundSpeciesDefinitionString."\n#undef _".$_."_BACKGROUND_SPEC_\n#define _".$_."_BACKGROUND_SPEC_  $cnt\n";
        $cnt++;
      }
      
      ampsConfigLib::AddLine2File($BackgroundSpeciesDefinitionString,"pic/picSpeciesMacro.dfn");
      
      #create the background Atmosphere Mass Table
      my @MassTable;
      
      foreach (@BackgroundSpeciesList) {
        push(@MassTable,"_".$_."__MASS_");
      }
      
      #create the conversion table between background and model speces (if such conversion is possible)
      my @Background2ModelSpeciesConversionTable;
      
      foreach (@BackgroundSpeciesList) {
        my $bspec=$_;
        my $nspec=-1;
        
        for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
          if ($bspec eq $SpeciesList[$i]) {
            $nspec=$i;
            last;
          }
        }
        
        push(@Background2ModelSpeciesConversionTable,$nspec);
      }
      
      #add new values of the variables into the code
      ampsConfigLib::ChangeValueOfArray("static const double BackgroundSpeciesMassTable\\[\\]",\@MassTable,"pic/pic.h");
      ampsConfigLib::ChangeValueOfVariable("static const int nTotalBackgroundSpecies",$nBackgroundSpecies,"pic/pic.h");
      ampsConfigLib::ChangeValueOfArray("static const int Background2ModelSpeciesConversionTable\\[\\]",\@Background2ModelSpeciesConversionTable,"pic/pic.h");
         
    }
    elsif ($s0 eq "COLLISIONMODE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_","pic/picGlobal.dfn");
      }
      else {
        die "Cannot recognize the option\n";
      }
    }
    elsif ($s0 eq "COLLISIONCROSSSECTION") {
      ($s0,$s1)=split(' ',$s1,2);
      
      my @CrossSection;
        
      for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
        my @tmp=((0)x$nBackgroundSpecies);
        push(@CrossSection, [@tmp]);
      }      
      
      if ($s0 eq "CONST") {
        my (@localModelSpeciesList,@localBackgroundSpeciesList,@localCrossectionList);
        my ($nBackgroungSpec,$nModelSpec);
                       
        $s1=~s/[()]/ /g;
        
        while (defined $s1) {
          ($s0,$s1)=split(' ',$s1,2);
          
          if ($s0 eq "CONST") {
            $nBackgroungSpec=0;
            $nModelSpec=-1;
            
            my ($t1,$t2,$t3,$cnt,$flag);
            ($t1,$t2,$t3,$s1)=split(' ',$s1,4);
            push(@localCrossectionList,$t3);
            
            $cnt=0;
            $flag=0;
            
            foreach (@BackgroundSpeciesList) {
              if ($_ eq $t2) {
                $nBackgroungSpec=$cnt;
                $flag=1;
                last;
              }
              
              $cnt++;
            }
            
            if ($flag eq 0) {
              die "Cannot recognize background specie $t2, const(ModelSpecies,BackgroundSpecies), the definition of the background species MUST be at the begining of the 'BackgroundSpeciesBlock'\n";
            }
            
            $cnt=0;
            $flag=0;
            
            foreach (@SpeciesList) {
              if ($_ eq $t1) {
                $nModelSpec=$cnt;
                $flag=1;
                last;
              }
              
              $cnt++;
            }
            
            if ($flag eq 0) {
              die "Cannot recognize model specie $t1, const(ModelSpecies,BackgroundSpecies), the definition of the background species MUST be at the begining of the 'BackgroundSpeciesBlock'\n";
            }
            
            $CrossSection[$nModelSpec][$nBackgroungSpec]=$t3;          
          }          
        }
        
        #print out the constant cross section table
        for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
          for (my $j=0;$j<$nBackgroundSpecies;$j++) {
            print "$CrossSection[$i]->[$j]\n";
          }
        }                
      }
      elsif ($s0 eq "FUNCTION") {
        my $FunctionName;
        
        $line=~s/=/ /g;
        ($FunctionName,$line)=split(' ',$line,2);
        ($FunctionName,$line)=split(' ',$line,2);
        
        while (defined $s1) {
          ($s0,$s1)=split(' ',$s1,2);
          ($FunctionName,$line)=split(' ',$line,2);
          
          if ($s0 eq "FUNCTIONNAME") {
            ($FunctionName,$line)=split(' ',$line,2);
            last;
          }            
        }
        
        #add the difinition of the function into the code
        ampsConfigLib::RedefineMacroFunction("_PIC_BACKGROUND_ATMOSPHERE__COLLISION_CROSS_SECTION_FUNCTION_","(t1,t2,t3,t4,t5,t6) ($FunctionName(t1,t2,t3,t4,t5,t6))","pic/pic.h");  
      }
      else {
        die "Cannot recognize the option\n";
      }
      
      #prepare the cross section definition array
      my ($t,@ConstantCrossSectionTable);
        
      for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
        $t="{";
          
        for (my $j=0;$j<$nBackgroundSpecies;$j++) {
          $t=$t."$CrossSection[$i][$j]";
            
          if ($j < ($nBackgroundSpecies - 1)) {$t=$t.",";}
        }
          
        push(@ConstantCrossSectionTable,$t."}");
      }
        
      #insert the constant cross section table into the code
      ampsConfigLib::ChangeValueOfArray("static const double BackgroundAtmosphereConstantCrossSectionTable\\[PIC::nTotalSpecies\\]\\[nTotalBackgroundSpecies\\]",\@ConstantCrossSectionTable,"pic/pic.h");
      
    }
           
    elsif ($s0 eq "COLLISIONSCATTERINGANGLE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ISOTROPIC") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_","_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_","pic/picGlobal.dfn");
      }
      elsif ($s0 eq "FUNCTION") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_","_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__USER_DEFINED_","pic/picGlobal.dfn");
      }
    }
  
    elsif ($s0 eq "INJECTCONDITIONBACKGROUNDPARTICLE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_PARTICLE_ACCEPTANCE_MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_PARTICLE_ACCEPTANCE_MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }
    }

    elsif ($s0 eq "REMOVECONDITIONMODELPARTICLE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__MODEL_PARTICLE_REMOVAL_MODE_","_PIC_MODE_ON_","pic/picGlobal.dfn");
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__MODEL_PARTICLE_REMOVAL_MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
      }
    }
    



    elsif ($s0 eq "LOADUSERDEFINITIONS") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        ampsConfigLib::RedefineMacro("_PIC_BACKGROUND_ATMOSPHERE__LOAD_USER_DEFINITION__MODE_","_PIC_MODE_ON_","pic/pic.h");
        
        my $FileName;
        $line=~s/=/ /g;
        ($FileName,$line)=split(' ',$line,2);
        ($FileName,$line)=split(' ',$line,2);        
        ($FileName,$line)=split(' ',$line,2);
        ($FileName,$line)=split(' ',$line,2);
        
        ampsConfigLib::RedefineMacro("_PIC_BACKGROUND_ATMOSPHERE__UDER_DEFINITION_","\"$FileName\"","pic/pic.h");
      }


    }

        
    elsif ($s0 eq "#ENDBACKGROUNDSPECIES") {
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
  }  
}

#=============================== Read UnomolecularReaction Block =============================
sub ReadUnimolecularReactions {
  
  #parameters of individual reaction
  my (@SourceSpecies,@ProductSpecies,@ReactionId,@LifeTime,@nProducts,$np,$id,$lt);
  
  #Global parameters   
  my $ReactionDescriptorString;
  my $UserLifetimeFunction;
  my $ReactionFlag=true;
  my $LifeTimeMultiplier=1;
  my $ReactionProcessor="";
  my $nMaxProducts=0;
  my $nTotalReactions=0;
  
  my ($s0,$s1,$s2,$spec);
    
  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[(=,)]/ /g;
    
    #get the first word in the sequence
    ($s0,$s1)=split(' ',$InputLine,2);
    
    
    if ($s0 eq "UNIMOLECULARREACTIONS") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        $ReactionFlag=true;
      }
      elsif ($s0 eq "OFF") {
        $ReactionFlag=false;
      }
      else {
        die "Cannot recognize the option (line=$InputLine, nline=$InputFileLineNumber)\n";
      }
    }
    elsif ($s0 eq "REACTIONLIFETIMEMULTIPLIER") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/\s+$//;
      $InputLine=~s/=/ /g;
      ($s0,$s1)=split(' ',$InputLine,2);
      ($LifeTimeMultiplier,$s1)=split(' ',$s1,2);
    }
    elsif ($s0 eq "PROCESSOR") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/\s+$//;
      $InputLine=~s/=/ /g;
      ($s0,$s1)=split(' ',$InputLine,2);
      ($ReactionProcessor,$s1)=split(' ',$s1,2);
    }
    elsif ($s0 eq "LIFETIMEUSERFUNCTION") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/\s+$//;
      $InputLine=~s/=/ /g;
      ($s0,$s1)=split(' ',$InputLine,2);
      ($UserLifetimeFunction,$s1)=split(' ',$s1,2);
    }
    
    elsif ($s0 eq "REACTION") {
      ($s0,$s1)=split(' ',$s1,2);
      
      $np=0;
      $id=-1;
      $nTotalReactions++;
      
      while (defined $s0) {
        if ($s0 eq "SOURCE") {
          ($s0,$s1)=split(' ',$s1,2);
          $spec=ampsConfigLib::GetElementNumber($s0,\@SpeciesList);
          push(@SourceSpecies,$spec);
        }
        elsif ($s0 eq "PRODUCT") {
          ($s0,$s1)=split(' ',$s1,2);
          $s0=~s/#/ /g;
          
          while (defined $s0) {
            ($s2,$s0)=split(' ',$s0,2);
            $spec=ampsConfigLib::GetElementNumber($s2,\@SpeciesList);
            
            if ($s2 eq "NONE") {
              #do nothing
            }
            elsif ($spec == -1) {
              die "Cannot recognize specie $s2 ($line) in $InputFileName.Assembled (Unimolecular Block)\n";
            }
            else {
              push(@ProductSpecies,$spec);
              $np++;
            }
          }
        }
        elsif ($s0 eq "ID") {
          ($id,$s1)=split(' ',$s1,2);
          push(@ReactionId,$id);
        }
        elsif ($s0 eq "LIFETIME") {
          ($lt,$s1)=split(' ',$s1,2);
          push(@LifeTime,$lt);
        }
        else {
          die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
        }
        
        ($s0,$s1)=split(' ',$s1,2);
      }
        
        #get the number of the reaction products
        push(@nProducts,$np);
        
        if ($np>$nMaxProducts) {
          $nMaxProducts=$np;
        }
        
    }
        

        
    elsif ($s0 eq "#ENDUNIMOLECULARREACTIONS") {
      my $n;
      my $cnt=0;     
      my @nSpeciesReactionNumber=((0)x$TotalSpeciesNumber);
      
      
      #output the line that describes the reactions
      print "$nMaxProducts\n";
      
      for ($n=0;$n<$nTotalReactions;$n++) {
        $nSpeciesReactionNumber[$SourceSpecies[$n]]++;
        
        if ($n>0) {
          $ReactionDescriptorString=$ReactionDescriptorString.",";
        }
        
        $ReactionDescriptorString=$ReactionDescriptorString."{".$ReactionId[$n].",".$LifeTime[$n].",".1.0/$LifeTime[$n].",".$SourceSpecies[$n].",".$nProducts[$n].",{";
        
        for ($spec=0;$spec<$nMaxProducts;$spec++) {
          if ($spec<$nProducts[$n]) {      
            $ReactionDescriptorString=$ReactionDescriptorString.$ProductSpecies[$cnt];
            $cnt++;
          }
          else {
            $ReactionDescriptorString=$ReactionDescriptorString."-1";
          }
          
          if ($spec<$nMaxProducts-1) {
            $ReactionDescriptorString=$ReactionDescriptorString.",";
          }
          else {
            $ReactionDescriptorString=$ReactionDescriptorString."}";
          }
        }
        
        $ReactionDescriptorString=$ReactionDescriptorString."}";
        
        print "$ReactionDescriptorString\n";
      }
      
      
      ampsConfigLib::ChangeValueOfVariable("static const int nTotalUnimolecularReactions",$nTotalReactions,"pic/pic.h"); 
      ampsConfigLib::ChangeValueOfVariable("static const int nMaxUnimolecularReactionProducts",$nMaxProducts,"pic/pic.h");
      ampsConfigLib::ChangeValueOfVariable("static const cUnimoleculecularReactionDescriptor UnimoleculecularReactionDescriptor\\[nTotalUnimolecularReactions\\]","{".$ReactionDescriptorString."}","pic/pic.h");
      
      #get the list of reaction each species is participating in
      my $nMaxSpeciesReactionNumber=0;
      my @TotalSpecieReactionRate=((0)x$TotalSpeciesNumber);
      my $SpeciesReactionList="";
#      my $cnt;
      
      for ($spec=0;$spec<$nTotalReactions;$spec++) {
        if ($nSpeciesReactionNumber[$spec]>$nMaxSpeciesReactionNumber) {
          $nMaxSpeciesReactionNumber=$nSpeciesReactionNumber[$spec];
        }
      }
      
      for ($spec=0;$spec<$nTotalReactions;$spec++) {
        if ($spec!=0) {
          $SpeciesReactionList=$SpeciesReactionList.",";
        }

        $SpeciesReactionList=$SpeciesReactionList."{";
        $cnt=0;
       
        for ($n=0;$n<$nTotalReactions;$n++) {
          if ($SourceSpecies[$n]==$spec) {
            if ($cnt!=0) {
              $SpeciesReactionList=$SpeciesReactionList.",";
            }
            
            $SpeciesReactionList=$SpeciesReactionList.$n;
            $TotalSpecieReactionRate[$spec]+=1.0/$LifeTime[$n];
            $cnt++;
          }
        }
        
        for ($n=$cnt;$n<$nMaxSpeciesReactionNumber;$n++) {
          if ($n!=0) {
            $SpeciesReactionList=$SpeciesReactionList.",";
          }
            
          $SpeciesReactionList=$SpeciesReactionList."-1";          
        }
      
        $SpeciesReactionList=$SpeciesReactionList."}";      
      }
      
      my $TotalReactionRateString="";
      
      for ($spec=0;$spec<$nTotalReactions;$spec++) {
        if ($spec!=0) {
          $TotalReactionRateString=$TotalReactionRateString.",";
        }
        
        $TotalReactionRateString=$TotalReactionRateString.$TotalSpecieReactionRate[$spec];
      }
      
      ampsConfigLib::ChangeValueOfVariable("static const int nMaxSpeciesUnimolecularReactionNumber",$nMaxSpeciesReactionNumber,"pic/pic.h");
      ampsConfigLib::ChangeValueOfVariable("static const int SpeciesUnimolecularReactionList\\[PIC::nTotalSpecies\\]\\[nMaxSpeciesUnimolecularReactionNumber\\]","{".$SpeciesReactionList."}","pic/pic.h");
      ampsConfigLib::ChangeValueOfVariable("static const double TotalSpecieUnimolecularReactionRate\\[PIC::nTotalSpecies\\]","{".$TotalReactionRateString."}","pic/pic.h");
      
      last;
    }
    else {      
      $line=~s/ //g;
      chomp($line);
   
      if (($line ne "") && (substr($line,0,1) ne '!')) {
        die "Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
      }
    }
  }  
  
}


#=============================== Read Species Block =============================
sub ReadSpeciesBlock {
  my $nspec=-1;
  
  my $s0;
  my $s1;
  
  my @MassArray=(0)x$TotalSpeciesNumber;
  my @ElectricChargeArray=(0)x$TotalSpeciesNumber;
  my @SpeciesFoundFlag=(0)x$TotalSpeciesNumber;
  
  my $SkipSpecieFlag=0;
  
  while ($line=<InputFile>) {
    ($InputFileLineNumber,$FileName)=split(' ',$line);
    $line=<InputFile>;
    
    ($InputLine,$InputComment)=split('!',$line,2);
    $InputLine=uc($InputLine);
    chomp($InputLine);
    $InputLine=~s/\s+$//; #remove spaces from the end of the line
 
    #substitute separators by 'spaces'
    $InputLine=~s/[=,]/ /g;
    
    #get the first word in the sequence
    ($s0,$s1)=split(' ',$InputLine,2);
    
    if ($s0 eq "#COMPONENT") {
      ($s0,$s1)=split(' ',$s1,2);
      
      $nspec=0;
      $SkipSpecieFlag=0;
      
      foreach (@SpeciesList) {
        if ($SpeciesList[$nspec] eq $s0) {
          $SkipSpecieFlag=1; #the species has been found in the list of those used in the simulation
          $SpeciesFoundFlag[$nspec]=1;
          last;
        }
        else {
          $nspec++;
        }      
      }
    }
    

    #pars the species properties    
    if ($SkipSpecieFlag == 1) {
      if ($s0 eq "MASS") {
        ($s0,$s1)=split(' ',$s1,2);
        $MassArray[$nspec]=$s0;
      }
      elsif ($s0 eq "CHARGE") {
        ($s0,$s1)=split(' ',$s1,2);
        $ElectricChargeArray[$nspec]=$s0."*ElectronCharge";
      }      
    }

    if ($s0 eq "#ENDSPECIES") {
      last;
    }
  }
  
  #add to the code the content of the Species Block
  
  
  #check is all modeled species where found
  for ($nspec=0;$nspec<$TotalSpeciesNumber;$nspec++) {
    if ($SpeciesFoundFlag[$nspec] == 0) {
      die "Cannot find physical parameters for $SpeciesList[$nspec]\n";
    }
  }

  
  #change appropriate variables in the code
  my @t;
  
  @t=@SpeciesList;
  foreach (@t) {
    $_="\"".$_."\"";
  }
  
  ampsConfigLib::ChangeValueOfArray("static const char ChemTable\\[\\]\\[_MAX_STRING_LENGTH_PIC_\\]",\@t,"pic/pic.h");
  ampsConfigLib::ChangeValueOfArray("static const double MolMass\\[\\]",\@MassArray,"pic/pic.h");
  ampsConfigLib::ChangeValueOfArray("static const double ElectricChargeTable\\[\\]",\@ElectricChargeArray,"pic/pic.h");
  
  #init the array of species type descriptors
  my @SpcecieTypeTable=("_PIC_SPECIE_TYPE__GAS_")x$TotalSpeciesNumber;
  ampsConfigLib::ChangeValueOfArray("static const int SpcecieTypeTable\\[\\]",\@SpcecieTypeTable,"pic/pic.h");

=comment
  #the chemical symbol table  
  print "const static char ChemicalSymbolsArray[$TotalSpeciesNumber]={\n";
  for ($nspec=0;$nspec<$TotalSpeciesNumber;$nspec++) {
    if ($nspec==$TotalSpeciesNumber-1) {
      print "\"$SpeciesList[$nspec]\"\n};\n";
    }
    else {
      print "\"$SpeciesList[$nspec]\",\n";
    }
  }
  
  #the mass table  
  print "const static double MassArray[$TotalSpeciesNumber]={\n";
  for ($nspec=0;$nspec<$TotalSpeciesNumber;$nspec++) {
    if ($nspec==$TotalSpeciesNumber-1) {
      print "\"$MassArray[$nspec]\"\n};\n";
    }
    else {
      print "\"$MassArray[$nspec]\",\n";
    }
  }
=cut


  #save the list of the species
  open(SPECLIST,">$InputFileName.Assembled.Species");
  my $size = scalar @SpeciesList;
  print SPECLIST "$size\n";
  
  foreach (@SpeciesList) {
    print SPECLIST "$_\n";
  }
  close (SPECLIST);
}



#=============================== Process AMPS' settings from .ampsConfig
sub ampsConfigSettings {
  if (-e ".ampsConfig.Settings") {
    my @Settings;
  
    open (AMPSSETTINGS,".ampsConfig.Settings") || die "Cannot open file\n";
    @Settings=<AMPSSETTINGS>;
    close (AMPSSETTINGS);
  
    
    foreach (@Settings){
      if (/^SPICEKERNELS=(.*)$/i) {ampsConfigLib::ChangeValueOfVariable("const char SPICE_Kernels_PATH\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$1."\"","models/exosphere/Exosphere.h"); next};
      if (/^SPICE=(.*)$/i) {`echo "SPICE=$1" >> Makefile.local`; next};
      
    }    
  }
}

=comment
#=============================== Change a value of a variable in the code  =============================
sub ChangeValueOfVariable {
  my $Variable=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    if ($_=~/$Variable/) {
      $Variable=~s/\\//g;
      $_=$Variable."=".$Value.";\n";
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
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro /) {
#      print "\#define $Macro $Value\n";
      $_="\#define $Macro $Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

sub RedefineMacroFunction {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro/) {
#      print "\#define $Macro$Value\n";
      $_="\#define $Macro$Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

sub AddMacro {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  open (FILEIN,">>$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";   
  print FILEIN "#undef $Macro\n#define $Macro $Value\n\n";
  close (FILEIN);
}
  
#=============================== Change a value of a array in the code  =============================
sub ChangeValueOfArray {
  my $Variable=$_[0];
  my @Value=@{$_[1]};
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
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

#===============================  Substitute a line in a source file =============================
sub SubstituteCodeLine {
  my $oldLinstKey=$_[0];
  my $newLine=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    if ($_=~m/($oldLinstKey)/) {
      $_=$newLine;
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

#===============================  RECURSIVE SUBSTITUDE OF A STRING IN THE SOURCE DIRECTORY ==========
sub RecursiveSubstitute {
  my $init=$_[0];
  my $final=$_[1];
  my $dir=$_[2];
  
  my $fname;
  my @FileList;
  
  opendir(DIR,$dir) or die "Cannot open directory $dir\n";
  @FileList=readdir(DIR);
  closedir(DIR);
  
  foreach $fname (@FileList) {
    #check if '$fname' is a directory
    
    if (($fname eq '.')||($fname eq '..')||($fname eq 'CVS')||($fname eq 'RCS')) {
      next;
    }
    
    if (-d "$dir/$fname") {
      #fname is a directory -> call the function recursevely
      RecursiveSubstitute($init,$final,"$dir/$fname");
    }
    else {
      #fname is a file -> make the substitution
      my @lines;
      
      open (FILE,"<$dir/$fname") || die "Cannot open file $dir/$fname\n";
      @lines=<FILE>;
      close (FILE);
      
      foreach (@lines) {
        $_=~s/$init/$final/g;
      }
      
      open (FILE,">$dir/$fname") || die "Cannot open file $dir/$fname\n";
      print FILE @lines;
      close (FILE);
    }
  }
}

=cut


