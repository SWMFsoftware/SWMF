#!/usr/bin/perl -i
use strict;

our @Arguments = ('src/Makefile.def', @ARGV);

my $config     = "share/Scripts/config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}
