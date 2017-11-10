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

my $source;
foreach $source (@source){

    # read source file into an array of lines
    open(FILE, $source);
    my @lines = <FILE>;
    close(FILE);

    # Replace old variable names with new names
    foreach (@lines){
	s/\bBLKtest\b/iBlockTest/gi;
	s/\bPROCtest\b/iProcTest/gi;
	s/\bVARtest\b/iVarTest/gi;
	s/\bDIMtest\b/iDimTest/gi;
	s/\biBLK\b/iBlock/gi;
	s/\boktest(_me)?\b/DoTest/gi;
    }

    my $text = join('', @lines);
    my $usetest = "  use BATL_lib, ONLY: &\n       test_start, test_stop";

    my $testvar;
    foreach $testvar (@testvar){
	$usetest .= ", $testvar" if $text =~ /\b$testvar\b/i;
    }

    my $subroutine;
    my $addnamesub;
    my $adddotest;
    my $testarguments;
    my $remove;
    foreach (@lines){

	# Add use statement for all test variables
	if(s/^Module/module/i){
	    $_ .= "\n$usetest\n";
	}

	# Fix separator line(s)
	if(/^(\s+\!)=+\s+$/){
	    $_  = $1 . "=" x (79 - length($1)) . "\n";
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
	$_ = "" if /parameter::\s*NameSub\s*=/ or /\s+logical\s?::\s?DoTest/
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
