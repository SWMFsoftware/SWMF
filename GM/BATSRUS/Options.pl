#!/usr/bin/perl -s

# Put switches into properly named variables
my $Help         = $h or $H or $help or $Help;
my $Equation     = $e or $equation;
my $UserModule   = $u or $user or $usermod;
my $Show         = $s or $show;

use strict;

# Print help if required
&print_help if $Help;

# Set error string for error messages
my $ERROR = 'Options.pl ERROR:';

# Define the main directories and files
my $Src='src';
my $SrcUser='srcUser';
my $UserMod="$Src/ModUser.f90";

my $EquationMod="$Src/ModVarIndexes.f90";

# Set or list the equations
&set_equation if $Equation;

# Set or list the user modules
&set_user_module if $UserModule;

# Show current options
&show_options if $Show;

#############################################################################

sub set_equation{

    if($Equation eq '1'){
	my @EquationModules;
	chdir $Src;
	@EquationModules = sort(glob("ModVarIndexes_*.f90"));
	for (@EquationModules){s/^ModVarIndexes_//; s/\.f90$//;}
	print "Available Equations:\n   ",join("\n   ",@EquationModules),"\n";
	chdir "..";
    }
}

#############################################################################

sub set_user_module{

    if($UserModule eq '1'){
	my @UserModules;
	chdir $SrcUser;
	@UserModules = ("Default", sort(glob("*.f90")));
	for (@UserModules){s/^ModUser//; s/\.f90$//;}
	print "Available User Modules:\n   ",join("\n   ",@UserModules),"\n";
	chdir "..";
	return;
    }

    my $File;
    if($UserModule eq "Default"){
	$File = "$Src/ModUserDefault.f90";
    }else{
	$File = "$SrcUser/ModUser$UserModule.f90";
    }
    die "$ERROR File $File does not exist!\n" unless -f $File;
    `cp $File $UserMod`; # die "$ERROR Could not cp $File $UserMod\n";
}

#############################################################################

sub show_options{

    open(FILE, $UserMod) or die "$ERROR Could not open $UserMod\n";
    my $Module='???';
    my $Version='???';
    while(<FILE>){
	if(/NameUserModule/){
	    $Module = $';
	    $Module = <FILE> if $Module =~ /\&\s*$/; # read continuation line
	    $Module =~ s/\s*=\s*//;   # remove equal sign
	    $Module =~ s/^\s*[\'\"]//;   # remove leading quotation mark
	    $Module =~ s/[\'\"].*\n//;  # remove trailing quotation marks
	}
	$Version = $1 if /VersionUserModule\s*=\s*([\d\.]+)/;
    }
    close(FILE);
    print "User Module = $Module, ver $Version\n";

    open(FILE, $EquationMod) or die "$ERROR Could not open $EquationMod\n";
    my $Equation='???';
    while(<FILE>){
	next unless /NameEquation\s*=\s*[\'\"]([^\'\"]*)/;
	$Equation = $1; last;
    }

    print "Equation    = $Equation\n";
}

#############################################################################

sub print_help{

    print "
Purpose:

    Set various options for BATSRUS, such as equations and user module.

Usage:

    Options.pl [-h] [-s] [-e[=EQUATION]] [-u[=USERMODULE]]

  -h              Print help message and stop

  -e              List all available equations.

  -e=EQUATION     Select equation EQUATION. 

  -u              List all the available user modules.

  -u=USERMODULE   Select the user module USERMODULE. 

  -s or -show     Show current settings.

Examples:

  Show current settings:

Options.pl -s

  List available options for equations and user modules:

Options.pl -e=? -u=?

  Select the MHD equations and the Default user module:

Options.pl -e=MHD -u=Default",
"\n\n";
    exit 0;
}
