#!/usr/bin/perl
#$Id$
#reads the enceladus--nultiplume model section of the input file

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);
use Class::Struct;

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

#define structures describing TigerStripes and individual plumes
struct cTigerStripe => {
  ID => '$',
  ActiveFlag => '$',
  TemparatureTable => '@',
  BulkVelocityTable => '@',
  SourceRateTable => '@'
};

#allocate a Tiger Stripe table and object
my @TigerStripeTable;
#%my $TigerStripe=cTigerStripe->new();

my %TigerStripe;
my $TigerStripeLine;
my $TigerStripeID;
my $nActiveTigerStripes=0;


$TigerStripe{ID}="";
$TigerStripe{ActiveFlag}=0;
$TigerStripe{TemperatureTable}=[((0)x$TotalSpeciesNumber)]; ##@UniformSourceRate]; #((0)x$TotalSpeciesNumber);
$TigerStripe{BulkVelocityTable}=[((0)x$TotalSpeciesNumber)];
$TigerStripe{SourceRateTable}=[((0)x$TotalSpeciesNumber)];


#$TigerStripe{TemparatureTable}[2]=314111;

#print"2:: $TigerStripe{TemparatureTable}[2]\n";

#$TigerStripe->ID("");
#$TigerStripe->ActiveFlag(0);
#$TigerStripe->{TemparatureTable}=((0)x$TotalSpeciesNumber);
#$TigerStripe->{BulkVelocityTable}=((0)x$TotalSpeciesNumber);
#$TigerStripe->{SourceRateTable}=((0)x$TotalSpeciesNumber);




my $t=1;

#$TigerStripe->TemparatureTable(0,2);

print "PASS INIT LINE\n";
#print "1: @{TigerStripe->TemparatureTable}->[0]\n";

while ($line=<InputFile>) {
  #re-init the Tiger Stripe and individual plume objects
  $TigerStripe{ActiveFlag}=0;
  
  for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
    $TigerStripe{TemparatureTable}[$i]=0;
    $TigerStripe{BulkVelocityTable}[$i]=0;
  }
  
  #process the new line of the input file
  ($InputFileLineNumber,$FileName)=split(' ',$line);
  $line=<InputFile>;
  
  $LineOriginal=$line;
  chomp($line);
  
  
  ($InputLine,$InputComment)=split('!',$line,2);
  $InputLine=uc($InputLine);
  chomp($InputLine);
 
  $InputLine=~s/[=():]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  
  
  
  #begin parcing of the input file
  #Tiger Stripe  
  if ($InputLine eq "TIGERSTRIPE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    
    
    while (defined $InputLine) {
      if ($InputLine eq "ID") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        $TigerStripe{ID}=$InputLine;
      }
      elsif ($InputLine eq "ACTIVE") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "ON") {
          $TigerStripe{ActiveFlag}=1;
        }
        elsif ($InputLine eq "OFF") {
          $TigerStripe{ActiveFlag}=0;
        }
      }
      elsif ($InputLine eq "TEMPERATURE") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "ALL") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          $TigerStripe{TemperatureTable}=[(($InputLine)x$TotalSpeciesNumber)];
        }
        else {
          my $id;
          
          $id=getSpeciesNumber($InputLine);
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          $TigerStripe{TemparatureTable}[$id]=$InputLine;
        }
      }
      elsif ($InputLine eq "BULKVELOCITY") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "ALL") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          $TigerStripe{BulkVelocityTable}=[(($InputLine)x$TotalSpeciesNumber)];
        }
        else {
          my $id;
          
          $id=getSpeciesNumber($InputLine);
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          $TigerStripe{BulkVelocityTable}[$id]=$InputLine;
        }
      }
      elsif ($InputLine eq "SOURCERATE") {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
        
        if ($InputLine eq "ALL") {
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          $TigerStripe{SourceRateTable}=[(($InputLine)x$TotalSpeciesNumber)];
        }
        else {
          my $id;
          
          $id=getSpeciesNumber($InputLine);
          ($InputLine,$InputComment)=split(' ',$InputComment,2);
          $TigerStripe{SourceRateTable}[$id]=$InputLine;
        }
      }  
      else {
        die "The option ($InputLine) is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
      }          
      
      if (defined $InputComment) {
        ($InputLine,$InputComment)=split(' ',$InputComment,2);
      }
      else {
        last;
      }      
    }
    
    #add the data entry to the Tiger Stripe line
    if ($nActiveTigerStripes!=0) {
      $TigerStripeLine=$TigerStripeLine.", ";
    }
    
      
    $TigerStripeID=$TigerStripeID."\n#undef _$TigerStripe{ID}__ID_\n#define _$TigerStripe{ID}__ID_ $nActiveTigerStripes\n";
    $nActiveTigerStripes++;
      
    $TigerStripeLine=$TigerStripeLine."{&TigerStripe_$TigerStripe{ID},{";
 
    #output the source rate
    for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
      if ($i!=0) {
        $TigerStripeLine=$TigerStripeLine.", ";
      }
      
      $TigerStripeLine=$TigerStripeLine."$TigerStripe{SourceRateTable}[$i]";
    }
    
    #output temeprature 
    $TigerStripeLine=$TigerStripeLine."},{";

    for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
      if ($i!=0) {
        $TigerStripeLine=$TigerStripeLine.", ";
      }
      
      $TigerStripeLine=$TigerStripeLine."$TigerStripe{TemperatureTable}[$i]";
    }

    #output bulk speed 
    $TigerStripeLine=$TigerStripeLine."},{";

    for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
      if ($i!=0) {
        $TigerStripeLine=$TigerStripeLine.", ";
      }
      
      $TigerStripeLine=$TigerStripeLine."$TigerStripe{BulkVelocityTable}[$i]";
    }


    #output ID and the Active Flag    
    $TigerStripeLine=$TigerStripeLine."},\"$TigerStripe{ID}\"";
    
    if ($TigerStripe{Active}==1) {
      $TigerStripeLine=$TigerStripeLine.",true}";
    }
    else {
      $TigerStripeLine=$TigerStripeLine.",false}";
    }
    
    
    print "LINE:::::\n $TigerStripeID\n";
    print "LINE:::::\n $TigerStripeLine\n";
    
  }


  #end statment of the block
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }
   
  else {
    die "Option ($InputLine) is unknown, line=$InputFileLineNumber ($InputFileName)\n";
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
	    die "Cannot find species $species\n";
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
