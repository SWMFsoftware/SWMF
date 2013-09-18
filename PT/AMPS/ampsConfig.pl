#!/usr/bin/perl
#$Id$
 
use strict;
use warnings;

use ampsConfigLib;

my $loadedFlag_MainBlock=0;
my $loadedFlag_SpeciesBlock=0;
my $loadedFlag_BackgroundSpeciesBlock=0;

my $InputFileName="moon.input";
my $line;
my $InputFileLineNumber=0;
my $InputLine;
my $InputComment;
my $s0;
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
  
}


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
  elsif ($InputLine eq "#SAMPLING") {
    Sampling();
  } 
  elsif ($InputLine eq "#USERDEFINITIONS") {
    UserDefinitions();
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
  
}

#modify the makefile for the compilation mode
my @makefilelines;
      
open (MAKEFILE,"<Makefile.def") || die "Cannot open Makefile.def\n";
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
      
open (MAKEFILEFILE,">Makefile.def") || die "Cannot open Makefile.def\n";
print MAKEFILEFILE @makefilelines;
close (MAKEFILEFILE);


#compile the code 
if ($CompileProcessedCodeFlag==1) {
  print "Compile the code\n";
  system("make");
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
      #substitute the following line in the Makefile.def 
      my $Variable;
      my $Argument;
      my @MakeFileContent;
      
      my $FoundVariable=0;
      
      ($line,$Argument)=split('!',$line,2);
      $Variable=$line;
      $Variable=~s/=/ /;
      ($Variable,$Argument)=split(' ',$Variable,2);
      ($Variable,$Argument)=split(' ',$Argument,2);
      
      if (-e "Makefile.def") {     
        open (MAKEFILE,"<","Makefile.def") || die "Cannot open Makefile.def\n";
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
          #the definition of the variable is not found -> add it to Makefile.def
          push(@MakeFileContent,$Variable."=".$Argument);
        }
      }
      
      #write the updated copy of the Makefile.def 
      open (MAKEFILE,">","Makefile.def");

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
  if (-e "Makefile.def") {
    my @MakeFileContent;
    my $FoundWSD=0;
    
    open (MAKEFILE,"<","Makefile.def") || die "Cannot open Makefile.def\n";
    @MakeFileContent=<MAKEFILE>;
    close(MAKEFILE);
    
    open (MAKEFILE,">","Makefile.def");

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
      open (MAKEFILE,"<","Makefile.def") || die "Cannot open Makefile.def\n";
      @MakeFileContent=<MAKEFILE>;
      close(MAKEFILE);
    
      open (MAKEFILE,">","Makefile.def");   
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
  
  #change prefix,error log file and diagnostic stream
  ampsConfigLib::RecursiveSubstitute('\$PREFIX:',$Prefix,$ampsConfigLib::WorkingSourceDirectory);
  ampsConfigLib::RecursiveSubstitute('\$ERRORLOG',$ErrorLog,$ampsConfigLib::WorkingSourceDirectory);
  ampsConfigLib::ChangeValueOfVariable("char PIC::DiagnospticMessageStreamName\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$DiagnosticStream."\"","pic/pic_init_const.cpp");
  ampsConfigLib::ChangeValueOfVariable("char PIC::OutputDataFileDirectory\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$OutputDirectory."\"","pic/pic_init_const.cpp");
  
  
  #update variables of the code distribution with the parameters of the "main block" of the input file 
  ampsConfigLib::ChangeValueOfVariable($InputFileVariable2CodeVariable{"TotalSpeciesNumber"},$TotalSpeciesNumber,$InputFileVariable2CodeSourceFile{"TotalSpeciesNumber"});
  
  #redefine the value of the macro that allows loading the user-defined table of the macros that describe the species used in the simulation 
  ampsConfigLib::RedefineMacro("_PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE_","_PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE__ON_","pic/picGlobal.dfn");

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
  
  #add the definitions of the species macros to 'UserDefinition.PIC.h'
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h\n";  
  @FileContent=<FILEIN>;
  close (FILEIN); 
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h\n";
  print FILEOUT "$MARKER__SPECIES_MACRO_DEFINIETION_USED_IN_SIMULATION\n";
  
  foreach (@FileContent) {
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
        
        
        ampsConfigLib::RedefineMacro("_PIC__PARTICLE_COLLISION_MODEL__MODE_","_PIC_MODE_OFF_","pic/picGlobal.dfn");
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
              
          if (($s0 eq -1) || ($s1 eq -1)) {
            die "Cannot recognize model specie (line=$InputFileLineNumber)\n";
          }
              
          $CrossSectionTable[$s0][$s1]=$t3;
          $CrossSectionTable[$s1][$s0]=$t3;                  
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
        
        if (-e "main/UserDefinition.PIC.h") {
          ampsConfigLib::RedefineMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__ON_","main/UserDefinition.PIC.h");
        }
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::RedefineMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_","pic/picGlobal.dfn");
        
        if (-e "main/UserDefinition.PIC.h") {
          ampsConfigLib::RedefineMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_","main/UserDefinition.PIC.h");
        }
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
      my @FileContent;
      my $cnt=0;
      
      open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h\n";  
      @FileContent=<FILEIN>;
      close (FILEIN); 
  
      open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/main/UserDefinition.PIC.h\n";
      
      foreach (@BackgroundSpeciesList) {
        print FILEOUT "#undef _"."$_"."_BACKGROUND_SPEC_\n";
        print FILEOUT "#define _"."$_"."_BACKGROUND_SPEC_  $cnt\n\n";
        $cnt++;
      }
      
      foreach (@FileContent) {
        print FILEOUT "$_";
      }
   
      close (FILEOUT);
      
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
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__ON_","main/UserDefinition.PIC.h");
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_MODE_","_PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_","main/UserDefinition.PIC.h");
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
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_","_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_","main/UserDefinition.PIC.h");
      }
      elsif ($s0 eq "FUNCTION") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_","_PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__USER_DEFINED_","main/UserDefinition.PIC.h");
      }
    }
  
    elsif ($s0 eq "INJECTCONDITIONBACKGROUNDPARTICLE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_PARTICLE_ACCEPTANCE_MODE_","_PIC_MODE_ON_","main/UserDefinition.PIC.h");
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_PARTICLE_ACCEPTANCE_MODE_","_PIC_MODE_OFF_","main/UserDefinition.PIC.h");
      }
    }

    elsif ($s0 eq "REMOVECONDITIONMODELPARTICLE") {
      ($s0,$s1)=split(' ',$s1,2);
      
      if ($s0 eq "ON") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__MODEL_PARTICLE_REMOVAL_MODE_","_PIC_MODE_ON_","main/UserDefinition.PIC.h");
      }
      elsif ($s0 eq "OFF") {
        ampsConfigLib::AddMacro("_PIC_BACKGROUND_ATMOSPHERE__MODEL_PARTICLE_REMOVAL_MODE_","_PIC_MODE_OFF_","main/UserDefinition.PIC.h");
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

#=============================== Read Species Block =============================
sub ReadSpeciesBlock {
  my $nspec=-1;
  
  my $s0;
  my $s1;
  
  my @MassArray=(0)x$TotalSpeciesNumber;
  
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
    }

 
    if ($s0 eq "#ENDSPECIES") {
      last;
    }
  }
  
  #add to the code the content of the Species Block
  
  
  #change appropriate variables in the code
  my @t;
  
  @t=@SpeciesList;
  foreach (@t) {
    $_="\"".$_."\"";
  }
  
  ampsConfigLib::ChangeValueOfArray("static const char ChemTable\\[\\]\\[_MAX_STRING_LENGTH_PIC_\\]",\@t,"pic/pic.h");
  ampsConfigLib::ChangeValueOfArray("static const double MolMass\\[nTotalSpecies\\]",\@MassArray,"pic/pic.h");
  
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


