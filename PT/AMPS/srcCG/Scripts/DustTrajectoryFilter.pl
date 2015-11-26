#!/usr/bin/perl

#select individual particle trajectories from AMPS' trajectory data file
#$Id$


#$ARGV[0]="amps.TrajectoryTracking.out=12.s=2.DUST:1.dat";

print "File name=$ARGV[0]\n";

open (FILE,"<$ARGV[0]") or die $!;
open (fSCRIPT,">script.mcr") or die $!;
open (fAcceptedTrajectories,">AcceptedTrajectories.dat") or die $!; 

print fSCRIPT "#!MC 1400\n"; 



#constants: direction of the axis, the radii and engle range
my $yAxis=0.0;
my $zAxis=1.0;
my $AngleRange=7.0;
my $DistanceOutOfPlane=0.1*2400.0E3;


my $rMin=0.00000001*2400.0E3;
my $rMax=3*2400.0E3;
my $AcceptanceRate=0.1; #0.001;

my $line;
my @s;
my $ZoneCounter=1;
my $ZoneONFlag=0;
my $nShowZones=0;
my $ZoneMustBeOFF=0;
my @ZoneData;

my $cosAngleRange=cos($AngleRange/180.0*3.141592654);

my $minCosAngle=cos(10.0/180.0*3.141592654);




$line=<FILE>; 
print "$line\n";
print fAcceptedTrajectories "$line\n";

$line=<FILE>;

#the direction at the beginig and the end of the trajectory
my @lInit=(0)x3;
my @x1=(0)x3;
my @x0=(0)x3; 
my $nPointsZone=0;

while ($line=<FILE>) {
   @s=split(' ',$line,5);

   if ($s[0] eq "ZONE") {
     #the previous zone is finished -> turn the zone off if needed

     #the zone sould be 'on' only when the direaction of the grains is reflected
     my $i;
     my $l1=0.0;
     my $l2=0.0;
     my $l3=0.0;
     
     for ($i=0;$i<3;$i++) {
       $l1+=$lInit[$i]**2;
       $l2+=($x1[$i]-$x0[$i])**2;
       $l3+=$lInit[$i]*($x1[$i]-$x0[$i]);
     }
          
     if (1==1) { #$nPointsZone>4) && ($lInit[0]/sqrt($l1)>cos(85.0/180.0*3.141592654)) && (-($x1[0]-$x0[0])/sqrt($l2) >0.0) ) {       #  ($l3/sqrt($l1*$l2)<$minCosAngle)) {
       
       if (rand()<$AcceptanceRate) {
         $ZoneONFlag=1;
       }
       else {    
         $ZoneONFlag=0;
       }
     }
     else {
       $ZoneONFlag=0;
     } 

     if (($ZoneONFlag == 1) && ($ZoneMustBeOFF==0) ) {
       print fSCRIPT "\$!ACTIVEFIELDMAPS += [$ZoneCounter]\n";
       $nShowZones+=1;
       
       print fAcceptedTrajectories "@ZoneData";
     }
     else {
       print fSCRIPT "\$!ACTIVEFIELDMAPS -= [$ZoneCounter]\n";
     } 

     #reset the flag
     $ZoneONFlag=0; 
     $ZoneMustBeOFF=0;
     $ZoneCounter+=1;
     $nPointsZone=0;
     
     undef(@ZoneData);
     push(@ZoneData,$line);
   }
   else {
     my $c;
     my $i;
     
     push(@ZoneData,$line);

if (sqrt($s[0]*$s[0]+$s[1]*$s[1]+$s[2]*$s[2])>20.0E3) { 
  
     if ($nPointsZone==0) {
       for ($i=0;$i<3;$i++) {
         $x0[$i]=$s[$i];
       }
     }
     elsif ($nPointsZone==1) {
       for ($i=0;$i<3;$i++) {
         $lInit[$i]=$s[$i]-$x0[$i];
         $x1[$i]=$s[$i];
       }
       
       if ($lInit[0]*$lInit[0]+$lInit[1]*$lInit[1]+$lInit[2]*$lInit[2]<1.0) {
         $nPointsZone=0;
       }
     } 
     else {
       my $length=0;

       for ($i=0;$i<3;$i++) {
         $length+=($s[$i]-$x1[$i])**2;
       }
       
       if (sqrt($s[0]*$s[0]+$s[1]*$s[1]+$s[2]*$s[2])>1010.0E3) {
         $ZoneMustBeOFF=1;
       }

       if ($length>1.0) {
         for ($i=0;$i<3;$i++) {
           $x0[$i]=$x1[$i];
           $x1[$i]=$s[$i];
         }
       }

     } 
     
     $nPointsZone++;
}
     

#=comment
#     if ($ZoneONFlag==0) {
#       my $r;
#
#       $r=sqrt($s[1]*$s[1]+$s[2]*$s[2]);
#
#       if ( (0.0<$s[0]) && ($s[0]<6000.0E3) && (abs($s[1])<$DistanceOutOfPlane) ) { ## ($rMin<$r) && ($r<$rMax) ) {
#
##       if ($r<5000.0E3) {
#
#         $c=abs($yAxis*$s[1]+$zAxis*$s[2])/$r;  
#
##        if ($c-$cosAngleRange>0) { 
#
#if (abs($s[1])<$DistanceOutOfPlane) {
#
#           if (rand()<$AcceptanceRate) {  
#             $ZoneONFlag=1;
#           }
#         }
#
#
#       }
#
#if (abs($s[1])>2000.0e3) {
#  $ZoneMustBeOFF=1;
#}
#
#     }
#
#if (abs($s[1])>2000.0e3) {
#  $ZoneMustBeOFF=1;
#}
#=cut

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
close (fAcceptedTrajectories);

#print statistics
print "Total number of zones: $ZoneCounter\n";
print "Zones ON: $nShowZones\n";


