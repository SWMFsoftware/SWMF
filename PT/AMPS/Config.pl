#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

our $Component       = 'PT';
our $Code            = 'amps';
our $MakefileDefOrig = 'Makefile.def.amps';
our @Arguments       = @ARGV;

my $Application;
my $MpiLocation;
my $Mpirun="mpirun";
my $TestRunProcessorNumber=4;

#create Makefile.local
if (! -e "Makefile.local") {
  `echo " " > Makefile.local`;
}


# Run the shared Config.pl script
my $config     = "share/Scripts/Config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

our %Remaining;   # Arguments not handled by share/Scripts/Config.pl

#sort the argument line to push the 'application' argument first
my $cntArgument=0;

foreach (@Arguments) {
  if (/^-application=(.*)/i) {
    if ($cntArgument!=0) {
      my $t;
      
      $t=$Arguments[0];
      $Arguments[0]=$Arguments[$cntArgument];
      $Arguments[$cntArgument]=$t;
      
      last;
    }
  }
  
  $cntArgument++;
}

#process the 'help' argument of the script
foreach (@Arguments) { 
   if ((/^-h$/)||(/^-help/)) { #print help message
     print "Config.pl can be used for installing and setting of AMPS\n";
     print "AMPS' settings can be defined in the SWMF's Config.pl: Config.pl -o=PT:spice-path=path,spice-kernels=path,ices-path=path,application=casename\n\n";
     print "Usage: Config.pl [-help] [-show] [-spice-path] [-spice-kernels] [-ices-path] [-application]\n";
     print "\nInformation:\n";
     print "-h -help\t\t\tshow help message.\n";
     print "-s -show\t\t\tshow current settings.\n";
     print "-spice-path=PATH\t\tpath to the location of the SPICE installation.\n";
     print "-spice-kernels=PATH\t\tpath to the location of the SPICE kernels.\n";
     print "-ices-path=PATH\t\t\tpath to the location of ICES\n";
     print "-application=case\t\tthe name of the model case to use.\n";
     
     exit;
   }
   
   
   if ((/^-s$/i)||(/^-show/)) { #print AMPS settings
     print "\nAMPS' Setting options:\n";
      
     if (-e ".ampsConfig.Settings") {
       my @f;
       
       open (FSETTINGS,"<.ampsConfig.Settings");
       @f=<FSETTINGS>;
       close(FSETTINGS);
       
       print @f;
       print "\n";
     }
     else {
       print "No custom settings\n";
     }
     
     exit;
   } 
}

#process the configuration settings
foreach (@Arguments) {
  if (/^-application=(.*)/i) {
    my $application = lc($1);
    
    `cp -f input/$application.* input/species.input .`;
    `echo "InputFileAMPS=$application.input" >> Makefile.local`;
      
    #copy all makefile settings from the input file into 'Makefile.local'
    #print "Config.pl modifies Makefile.local settings for application $1\n";
    open (INPUTFILE,"<","$application.input") || die "Cannot open $application.input\n";
      
    while (my $line=<INPUTFILE>) {
      my ($p0,$p1);
        
      chomp($line);
      ($p0,$p1)=split(' ',$line,2);
      $p0=lc($p0);
      $p0=~s/^\s+//;
        
      if ($p0 =~ /^makefile/) {
        `echo "$p1" >> Makefile.local`;
      }
    }
      
    close (INPUTFILE);   
    
    `echo "APPLICATION=$1" >> .ampsConfig.Settings`; 
    next
  };
    
  if (/^-mpi=(.*)$/i)        {$MpiLocation=$1;                next}; 
  if (/^-np=(.*)$/i)         {$TestRunProcessorNumber=$1;     next};
  if (/^-spice-path=(.*)$/i)      {`echo "SPICE=$1" >> Makefile.local`;`echo "SPICE=$1" >> .ampsConfig.Settings`;     next}; 
  if (/^-spice-kernels=(.*)$/i)    {`echo "SPICEKERNELS=$1" >> .ampsConfig.Settings`;     next};
  if (/^-ices-path=(.*)$/i)       {`echo "ICESLOCATION=$1" >> .ampsConfig.Settings`;     next};
  
  warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


exit 0;
