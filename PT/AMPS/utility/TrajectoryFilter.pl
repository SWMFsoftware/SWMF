#!/usr/bin/perl

print "File name=$ARGV[0]\n";

open (FILE,"<$ARGV[0]") or die $!;
open (fSCRIPT,">script.mcr") or die $!;
 

print fSCRIPT "#!MC 1400\n"; 


#constants: direction of the axis, the radii and engle range
my $yAxis=0.0;
my $zAxis=1.0;
my $AngleRange=7.0;
my $DistanceOutOfPlane=0.1*2400.0E3;


my $rMin=0.00000001*2400.0E3;
my $rMax=3*2400.0E3;
my $AcceptanceRate=1.0; #0.001;

my $line;
my @s;
my $ZoneCounter=1;
my $ZoneONFlag=0;
my $nShowZones=0;
my $ZoneMustBeOFF=0;

my $cosAngleRange=cos($AngleRange/180.0*3.141592654);


$line=<FILE>; 
print "$line\n";

$line=<FILE>;

while ($line=<FILE>) {
   @s=split(' ',$line,4);

   if ($s[0] eq "ZONE") {
     #the previous zone is finished -> turn the zone off if needed

     if (($ZoneONFlag == 1) && ($ZoneMustBeOFF==0) ) {
       print fSCRIPT "\$!ACTIVEFIELDMAPS += [$ZoneCounter]\n";
       $nShowZones+=1;
     }
     else {
       print fSCRIPT "\$!ACTIVEFIELDMAPS -= [$ZoneCounter]\n";
     } 

     #reset the flag
     $ZoneONFlag=0; 
     $ZoneMustBeOFF=0;
     $ZoneCounter+=1;
   }
   else {
     my $c;

     if ($ZoneONFlag==0) {
       my $r;

       $r=sqrt($s[1]*$s[1]+$s[2]*$s[2]);

       if ( (0.0<$s[0]) && ($s[0]<6000.0E3) && (abs($s[1])<$DistanceOutOfPlane) ) { ## ($rMin<$r) && ($r<$rMax) ) {

#       if ($r<5000.0E3) {

         $c=abs($yAxis*$s[1]+$zAxis*$s[2])/$r;  

#        if ($c-$cosAngleRange>0) { 

if (abs($s[1])<$DistanceOutOfPlane) {

           if (rand()<$AcceptanceRate) {  
             $ZoneONFlag=1;
           }
         }


       }

if (abs($s[1])>2000.0e3) {
  $ZoneMustBeOFF=1;
}

     }

if (abs($s[1])>2000.0e3) {
  $ZoneMustBeOFF=1;
}

   } 
}

#the last zone 
if ($ZoneONFlag == 0) {
  print fSCRIPT "\$!ACTIVEFIELDMAPS -= [$ZoneCounter]\n";
}
else {
  print fSCRIPT "\$!ACTIVEFIELDMAPS += [$ZoneCounter]\n";
}

#close the opne files 
close (FILE);
close (fSCRIPT);

#print statistics
print "Total number of zones: $ZoneCounter\n";
print "Zones ON: $nShowZones\n";


