#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#^CFG FILE CONFIGURE
# Find and run share/Scripts/Configure.pl -c=CFG
use strict;

my $ConfigureScript = "share/Scripts/Configure.pl";
my @ConfigureScript = ($ConfigureScript, 
		       "../../$ConfigureScript",
		       "../../../$ConfigureScript");

my $Script;
foreach $Script (@ConfigureScript){
    if(-f $Script){
	my $Command = "$Script -c=CFG ".join(' ',@ARGV);
	exec($Command) or die "Could not execute $Command\n";
    }
}

die "Could not find $ConfigureScript\n";   
