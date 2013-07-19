#!/usr/bin/perl -i
use strict;

my $Install;
my $Uninstall;
my $Compiler;
my $IsComponent;

foreach (@ARGV){
    if(/^-install(=.*)?$/)    {my $value=$1;
                               $IsComponent=1 if $value =~ /^=c/i;
                               $IsComponent=0 if $value =~ /^=s/i;
                               $Install=1;                     next;}
    if(/^-uninstall$/i)       {$Uninstall=1;                   next};
    if(/^-compiler=(.*)$/i)   {$Compiler=$1;                  next};
}

if($Install){

    $Compiler = "openmpicxx" unless $Compiler;

    if($IsComponent){
	`echo "include ../../Makefile.def" > Makefile.def`;
	`echo "include ../../Makefile.conf" > Makefile.conf`;
    }else{
	`cp -f Makefile.def.amps Makefile.def`;
	`echo "COMPILE.c = $Compiler" > Makefile.conf`;
    }
}

if($Uninstall){
    `make allclean`;
    `rm -f Makefile.def Makefile.conf`;
}

exit 0;
