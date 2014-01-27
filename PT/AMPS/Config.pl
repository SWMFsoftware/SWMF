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
        
      if ($p0 =~ /makefile/) {
        `echo "$p1" >> Makefile.local`;
      }
    }
      
    close (INPUTFILE);     
    next
  };
    
  if (/^-mpi=(.*)$/i)        {$MpiLocation=$1;                next}; 
  if (/^-np=(.*)$/i)         {$TestRunProcessorNumber=$1;     next};
  if (/^-spice=(.*)$/i)      {`echo "SPICE=$1" >> Makefile.local`;     next}; 
  if (/^-kernels=(.*)$/i)    {`echo "SPICEKERNELS=$1" >> .ampsConfig.Settings`;     next};
  if (/^-ices=(.*)$/i)       {`echo "ICESLOCATION=$1" >> .ampsConfig.Settings`;     next};
  
  warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


exit 0;
