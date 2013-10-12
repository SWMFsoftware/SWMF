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

# Run the shared Config.pl script
my $config     = "share/Scripts/Config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

our %Remaining;   # Arguments not handled by share/Scripts/Config.pl

foreach (@Arguments){
    if(/^-application=(.*)/i) {$Application=$1;                next};
    if(/^-mpi=(.*)$/i)        {$MpiLocation=$1;                next}; 
    if(/^-np=(.*)$/i)         {$TestRunProcessorNumber=$1;     next}; 

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

if($Application){
    my $application = lc($Application);
#   `rm -rf main; cp -r src$Application main; rm -rf main/CVS`;
    `cp -f input/$application.* input/species.input .`;
}

exit 0;
