#!/usr/bin/perl
#$Id$
#script for converting of the RMOC shape files into .mail and TECPLOT formats

my $ShapeFileName_RMOC="";

#read the argument line 
foreach (@ARGV) {
  if ((/^-h$/)||(/^-help/)) {
    print "The script converts RMOC shape files into two files that can be viewed with TECPLOT and loaded with AMPS\n";
    print "Argument list: -h -help -shape=[] is the name of the shape file\n";
    exit; 
  }
  elsif (/^-shape=(.*)/i) {
    $ShapeFileName_RMOC=$1;
  }
  else {
    die "Something is wrong with the argument list\n";
  }
} 

print "Name of the RMOC shape file: $ShapeFileName_RMOC\n";
open (fShape,"<$ShapeFileName_RMOC") || die "Cannot open file $ShapeFileName_RMOC\n";
open (fOutTECPLOT,">$ShapeFileName_RMOC.TECPLOT.dat");
open (fOutSHAPE,">$ShapeFileName_RMOC.mail");

#determine the number of the surface elements and the nodes 
my ($line,$s0,$s1,$s2,$s3,$s4,$i);
my $nTotalNodes=0;
my $nTotalFaces=0;

while ($line=<fShape>) {
  chomp ($line);
  $line=~s/=/ /g; 
  ($s0,$line)=split(' ',$line,2);
  
  if ($s0 eq "NUMBER_OF_VERTICES") {
    $nTotalNodes=$line;
  }
  elsif ($s0 eq "NUMBER_OF_FACES") {
    $nTotalFaces=$line;
    last;
  }
}

print "Total number of nodes: $nTotalNodes \nTotal Number of the surface elements: $nTotalFaces\n";

print fOutTECPLOT "VARIABLES=\"X\",\"Y\",\"Z\"\n"; 
print fOutTECPLOT "ZONE N=$nTotalNodes, E=$nTotalFaces, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n"; 

print fOutSHAPE "$nTotalNodes $nTotalFaces\n"; 

#skip two lines in the file 
$line=<fShape>;
$line=<fShape>;

#read the nodes of the surace mesh
for ($i=0;$i<$nTotalNodes;$i++) {
  $line=<fShape>;
  chomp ($line);
  ($s0,$s1,$s2,$s3)=split(' ',$line,4);

  $s1=1000.0*$s1;
  $s2=1000.0*$s2;
  $s3=1000.0*$s3;

  print fOutTECPLOT "$s1 $s2 $s3\n";
  print fOutSHAPE "$s1 $s2 $s3\n";
}
  
#read faces of the suraface mesh
for ($i=0;$i<$nTotalFaces;$i++) {
  $line=<fShape>;
  chomp ($line);
  ($s0,$s1,$s2,$s3,$s4)=split(' ',$line,4); 

  print fOutTECPLOT "$s2 $s3 $s4\n";
  print fOutSHAPE "$s2 $s3 $s4\n";
}

close (fOutTECPLOT);
close (fOutSHAPE);
