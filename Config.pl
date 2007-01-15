#!/usr/bin/perl -i
use strict;

our @Arguments  = ('SWMF','CON/Makefile.def',@ARGV);
my $config     = "share/Scripts/config.pl";
require $config or die "Could not find $config!\n";
