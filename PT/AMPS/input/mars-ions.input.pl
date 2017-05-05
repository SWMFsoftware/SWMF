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
    
  #the number, locations, energy range, and the number of the ebergy intervals used in the spherical sampling surfaces 
  if ($InputLine eq "SPHERICALSHELLS") {
    ($s0,$InputComment)=split(' ',$InputComment,2);
    
    if ($s0 eq "ON")  {
      my @ShellRadiiTable;
      my $ShellRadiiTableLength=0;
      
      #sampling will occurs
      ampsConfigLib::ChangeValueOfVariable("bool MarsIon::Sampling::SphericalShells::SamplingMode","true","main/mars-ions.cpp");
      
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
          ampsConfigLib::ChangeValueOfArray("double MarsIon::Sampling::SphericalShells::SampleSphereRadii\\[MarsIon::Sampling::SphericalShells::nSphericalShells\\]",\@ShellRadiiTable,"main/mars-ions.cpp");
          ampsConfigLib::ChangeValueOfVariable("const int nSphericalShells",$ShellRadiiTableLength,"main/mars-ions.h");                   
        }
        
        #check whether the entry is another setting parameter
        if ($s0 eq "EMIN") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("double MarsIon::Sampling::SphericalShells::minSampledEnergy",$s0,"main/mars-ions.cpp");         
        }
        elsif ($s0 eq "EMAX") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("double MarsIon::Sampling::SphericalShells::maxSampledEnergy",$s0,"main/mars-ions.cpp");          
        }
        elsif ($s0 eq "NLEVELS") {
          ($s0,$InputComment)=split(' ',$InputComment,2);
          ampsConfigLib::ChangeValueOfVariable("int MarsIon::Sampling::SphericalShells::nSampledLevels",$s0,"main/mars-ions.cpp");
        }
        else {
          die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
        }       
      }      
    }
    elsif ($s0 eq "OFF") {
       ampsConfigLib::ChangeValueOfVariable("bool MarsIon::Sampling::SphericalShells::SamplingMode","false","main/mars-ions.cpp");
    }     
    else {
      die "Cannot recognize $s0, line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
    }  
  }

  
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
  
}
