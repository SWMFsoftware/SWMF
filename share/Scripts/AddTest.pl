#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help = $h or $help;
my $Debug = $D;

# Allow in-place editing
$^I = "";

use strict;

my @source = @ARGV;

if($Help or not @source){
    print '
This script can be used to add calls of test_start and test_stop
subroutines from BATL_test and improve the formatting of the code.

Usage:

      AddTest.pl [-h] FILE1 [FILE2 ...]

-h    This help message

FILE1 FILE2 ... 
      Source files to be checked. Note that the script should be run in 
      the directory where the source files are.

Examples:

See the unused variables for all Fortran files:

     AddTest.pl *.f90

';
    exit;
}

# New test variable names
my @testvar = 
    ("iTest", "jTest", "kTest", "iBlockTest", "iProcTest", "iVarTest", "iDimTest", "xTest", "yTest", "zTest");

# Simple Fortran types with possible (len=..) and (kind=..) attributes:
my $SimpleType = '(real|integer|logical|character)(\s*\([^\)]+\))?\b';

# Obsolete Fortran types with *NUMBER, e.g. real*8 character*10
my $ObsoleteType = '(real|integer|logical|character)\s*\*\s*\d+';

# Derived Fortran type
my $DerivedType = 'type\s*\(\s*\w+\s*\)';

# Any Fortran Type
my $AnyType = "($SimpleType|$ObsoleteType|$DerivedType)";

my $indent;      # indentation at the beginning of the "unit"
my $unittype;    # program, subroutine or function
my $unitname;    # name of program unit
my $declaration; # true in declaration part
my $source;
foreach $source (@source){

    # read source file into an array of lines
    open(FILE, $source);
    my @lines = <FILE>;
    close(FILE);

    my $i = -1;
    foreach (@lines){
	# line index
	$i++;

	# Fix copyright message
	s/University of Michigan, portions used/University of Michigan,\n!  portions used/;

	# Replace old variable names with new names
	s/\bBLKtest\b/iBlockTest/gi;
	s/\bPROCtest\b/iProcTest/gi;
	s/\bVARtest\b/iVarTest/gi;
	s/\bDIMtest\b/iDimTest/gi;
	s/\biBLK\b/iBlock/gi;
	s/\boktest(_me)?\b/DoTest/gi;

	# remove original separator lines !==== and !-----
	$_ = '' if /^\s+\![\=\-]+\!?\s+$/;

	# Find start of subroutines and functions
	if(/^(\s*)(subroutine|function|$AnyType\s+function)\s(\w+)/i){
	    $indent = $1;
	    $unittype = $2;
	    $unitname = $3;
	    $declaration = 1;
	    next;
	}
	if($declaration){
	    # remove unnecessary "implicit none" from subroutines/functions
	    $_ = "" if /^\s+implicit\s+none/i;

	    # Skip empty lines and comments
	    next if s/^\s*$// or /^\s*\!/;

            # skip continuation lines
            next if $lines[$i-1] =~ /\&\s*(\!.*)?$/;
	    
	    # skip use statements
	    next if s/(\s*)use\b/$1use/i;

	    # skip variable declarations
	    next if /^\s*($AnyType)/i;

	    # End of declarations: insert !----- separator line
	    #print "end of declaration at line=$i for unit=$unitname, indent=$indent, line=$_";

	    $_ = $indent . '  !' . "-" x (76-length($indent)) . "\n" . $_;
	    $declaration = 0;
	}else{
	    s/\bcycle\b/CYCLE/;
	    s/\bexit\b/EXIT/;
	    s/\breturn\b/RETURN/;
	}

	# put in separator line at the end of methods
	if(/^(\s*)(contains|end\s+subroutine|end\s+function)\b/){
	    $_ .= $1 . "!" . ("=" x (78-length($1))) . "\n";
	    $unitname = "";
	    
	}

    }

    my $text = join('', @lines);
    my $usetest = "  use BATL_lib, ONLY: &\n       test_start, test_stop";

    my $testvar;
    foreach $testvar (@testvar){
	$usetest .= ", $testvar" if $text =~ /\b$testvar\b/i;
    }

    my $module;
    my $subroutine;
    my $addnamesub;
    my $adddotest;
    my $testarguments;
    my $remove;
    foreach (@lines){

	# Add use statement for all test variables at the beginning of the module
	if(s/^Module\s*(\w+)/module $1/i){
	    $module = $1;
	    $_ .= "\n$usetest\n";
	    print "Working on module $module\n";
	}

	if(/^  subroutine\s+(\w+)(.*)/){
	    $subroutine = $1;

	    if($2 =~ /\biBlock\b/i){
		$testarguments = "NameSub, DoTest, iBlock";
	    }else{
		$testarguments = "NameSub, DoTest";
	    }
	    next;
	}

	next unless $subroutine;

	# Remove original declarations for sake of uniformity
	$_ = "" if /parameter\s*::\s*NameSub\s*=/ or /\s+logical\s?::\s?DoTest/
	    or /^\s+call set_oktest/;

	# Remove "if(iProc==PROCtest .and. iBlock==BLKtest)...endif"
	$remove = 1
	    if /^\s+if.*then\s*$/ and /\biProcTest\b/ and /\biBlockTest\b/;

	if($remove){
	    # remove code
	    $remove = 0 if /\s+end\s?if/i;
	    $_ = "";
	    next;
	}

	# Declare DoTest, set NameSub, and add call test_start 
	# after the declarations
	if(/  !--------------------------------/){
	    chop;
	    $_ = "
    logical:: DoTest
    character(len=*), parameter:: NameSub = '$subroutine'
    !-------------------------------------------------------------------------
    call test_start($testarguments)

";
	}

	# Add call test_stop to the end of the subroutine
	if(/^  (contains|end subroutine)/){
	    $_ = "    call test_stop($testarguments)\n\n" . $_;
	    $subroutine = '';
	}
	
    }

    my $orig = $source."_orig_";
    `mv $source $orig` unless -f $orig;
    open(FILE, ">$source");
    print FILE @lines;
    close FILE;

}
