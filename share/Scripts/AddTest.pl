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

my $addnamesub;    # true if namesub should be declared
my $addtest;       # true if call test_start/stop should be added
my $namesubline;   # declaration of NameSub = ...
my $dotestline;    # declaration of DoTest logical
my $separatorline; # !---- line
my $indent;        # indentation at the beginning of the "unit"
my $unittype;      # program, subroutine or function
my $unitname;      # name of program unit
my $declaration;   # true in declaration part
my $usemodmain;    # true inside use ModMain... declaration
my $testarguments; # arguments for test_start and test_stop
my $removeoktest;  # true while removing old oktest code

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

	# remove trailing spaces
	s/\s+\n/\n/;

	# Fix copyright message
	s/Michigan, portions used/Michigan,\n!  portions used/;

	# Replace old variable names with new names
	s/\bBLKtest\b/iBlockTest/gi;
	s/\bPROCtest\b/iProcTest/gi;
	s/\bVARtest\b/iVarTest/gi;
	s/\bDIMtest\b/iDimTest/gi;
	s/\biBLK\b/iBlock/gi;
	s/\bnBLK\b/MaxBlock/gi;
	s/\boktest\b/DoTest/gi;
	s/\boktest_me\b/DoTestMe/gi;

	# remove original separator lines !==== and !-----
	$_ = '' if /^\s+\![\=\-]+\!?$/;

	# fix comments !blabla --> ! blabla
	s/^(\s*\![\!\$]+)([^\s\/\\\-\=\!\$])/$1 $2/;

	# Find start of subroutines and functions
	if(/^(\s*)(subroutine|function|$AnyType\s+function)\s+(\w+)/i){
	    $indent = $1;
	    $unittype = $2;
	    $unitname = $+;
	    $declaration = 1;

	    if($unittype eq "subroutine" and $indent eq "  "){
		$addtest = 1;
	    }else{
		$addtest = 0;
	    }

	    if(/\biBlock\b/i){
		$testarguments = "NameSub, DoTest, iBlock";
	    }else{
		$testarguments = "NameSub, DoTest";
	    }

	    $dotestline = "$indent  logical:: DoTest\n";
	    $namesubline = "$indent  ".
		"character(len=*), parameter:: NameSub = '$unitname'\n";
	    $separatorline = "$indent  !" . "-" x (76-length($indent)) . "\n";

	    #print "indent=$indent, unittype=$unittype, unitname=$unitname"
	    #	.", addtest=$addtest\n";

	    next;
	}
	if($declaration){
	    # remove unnecessary "implicit none" from subroutines/functions
	    $_ = "" if /^\s+implicit\s+none/i;

	    # Skip empty lines and comments
	    next if /^$/ or s/^\s*\n/\n/ or /^\s*\!/;

            # Fix ModSomething  ,only : --> ModSomething, ONLY:
            s/\bonly\s*:\b/ONLY:/i;
	    s/,ONLY/, ONLY/;
	    s/\s+(,\s+ONLY)/$1/;

            # remove iTest ... iVarTest from use ModMain
	    $usemodmain=1 if /^\s+use ModMain/i;
	    if($usemodmain){
		$usemodmain = 0 unless /\&$/;
		s/\b(i|j|k|iBlock|iProc|iVar|iDim|x|y|z)Test\b\s*,?//g;
		s/(,\s*)?(i|j|k|iBlock|iProc|iVar|iDim|x|y|z)Test$//g;
		next if s/^\s+use ModMain\s*,\s+ONLY:\s*$//i;
	    }

            # skip continuation lines
            next if $lines[$i-1] =~ /\&\s*(\!.*)?$/;

	    # skip use statements
	    next if s/(\s*)use\b/$1use/i;

	    # Remove old declarations of DoTest and DoTestMe
	    if(/^\s+logical\s*::\s*DoTest/){
		s/DoTest(Me)?\s*(=\s*\.(true|false)\.\s*)?,?\s*//g;
		$_ = '' if /^\s+logical\s*::\s*$/;
		$addtest = 1;
		next;
	    }

	    # Remove original NameSub declarations for sake of uniformity
	    if(/^\s+character.*parameter\s*::\s*NameSub\s*=/){
		$addnamesub = 1;
		$_ = "";
		next;
	    }

	    # skip all other variable declarations
	    next if /^\s*($AnyType)/i;

	    # End of declarations 
	    $declaration = 0;

	    # Declare DoTest if needed
	    $lines[$i-1] .= $dotestline if $addtest;

	    # Add NameSub declaration if needed
	    $lines[$i-1] .= $namesubline if $addnamesub or $addtest;
	    $addnamesub = 0;

	    # insert !----- separator line to the end of the next line
	    $lines[$i-1] .= $separatorline;

	    # add call test line
	    $lines[$i-1] .= $indent."  call test_start($testarguments)\n"
		if $addtest and not /\btest_start/;;
	}
	# fix code after the declaration part

	# Remove "if(iProc==iProcTest .and. iBlock==iBlockTest)...endif"
	$removeoktest = 1
	    if /^\s+if.*then$/ and /\biProcTest\b/ and /\biBlockTest\b/;

	if($removeoktest){
	    $removeoktest = 0 if /^\s+end\s?if/i;
	    $_ = "";
	    next;
	}
	    
	# Remove simple call set_oktest
	$_ = '' if /^\s+call\s+set_oktest/i;

	# DoTestMe --> DoTest
	s/\bDoTestMe\b/DoTest/gi;

	# Fix some common coding issues

	# Capitalize jumps
	s/\bcycle\b/CYCLE/;
	s/\bexit\b/EXIT/;
	s/\breturn\b/RETURN/;
	s/\bgo\s*to\b/GOTO/;

	# Obsolete relation operators
	s/\s*\.eq\.\s*/ == /ig;
	s/\s*\.ne\.\s*/ \/= /ig;
	s/\s*\.ge\.\s*/ \>= /ig;
	s/\s*\.le\.\s*/ \<= /ig;
	s/\s*\.gt\.\s*/ \> /ig;
	s/\s*\.lt\.\s*/ \< /ig;

	# Obsolete named constants
	s/\bcZero\b/0.0/ig;
	s/\bcHalf\b/0.5/ig;
	s/\bcOne\b/1.0/ig;

	# put in call test_stop() and separator line at the end of methods
	if(/^(\s*)(contains|end\s+subroutine|end\s+function)\b/){
	    $indent = $1;
	    $_ = $indent."  call test_stop($testarguments)\n".$_
		if $addtest and $lines[$i-1] !~ /\btest_stop/;

	    # extra indentation for separator after "contains" line
	    $indent .= "  " if /contains/;
	    $_ .= "$indent!" . ("=" x (78-length($indent))) . "\n";

	    $indent = "";
	    $unitname = "";
	    $addtest = 0;
	    $testarguments = "";
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
    my $i = -1;
    foreach (@lines){
	$i++;

	# Add use statement for all test variables at the beginning of module
	if(s/^Module\s*(\w+)/module $1/i){
	    $module = $1;
	    $_ .= "\n$usetest\n" unless $lines[$i+3] =~ /test_start/;
	    last;
	}
    }

    my $orig = $source."_orig_";
    `mv $source $orig` unless -f $orig;
    open(FILE, ">$source");
    print FILE @lines;
    close FILE;

}
