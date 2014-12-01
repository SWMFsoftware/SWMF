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


my $InputFileName="mars.input.Assembled.Block";   #$ARGV[0];
my $SpeciesFileName="mars.input.Assembled.Species";
my $WorkingSourceDirectory="srcTemp";  #$ARGV[1];
my $InputDirectory=$ARGV[2];  #the location of the data files that will be used in the model run

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

my $DataDirectory="data/input/Mars/";
my $DataCaseFiles=" ";

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


#add the model header to the pic.h
open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
print PIC_H "#include \"Mars.h\"\n";
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
  
  if ($InputLine eq "CASE") {
    my ($s0,$s1);
    
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=()]/ /g;
    ($s0,$s1,$InputLine)=split(' ',$InputLine,3);
    
    $s1=~s/ //g;
    $DataCaseFiles=$DataDirectory.$s1."/*.h";

    `mkdir -p $InputDirectory`;

    if ( -d $DataDirectory.$s1 ) {
	    `cp $DataCaseFiles $InputDirectory`;
    }
    else {
	     die "The case selected entitled $InputLine does not exist, line=$InputFileLineNumber ($InputFileName)\n";
    }
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
