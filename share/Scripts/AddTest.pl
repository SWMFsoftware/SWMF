#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $defaultexception = 
    "set_parameters|correct_electronfluid_efield_cell|select_fluid|".
    "set_yzr|user_get_b0|user_interface.f90|ModUserEmpty.f90";

my $Help       = ($h or $help);
my $Verbose    = ($v or $verbose);
my $Redo       = ($r or $redo);
my $Undo       = ($u or $undo);
my $TestLevel  = ($l or $level or 1);
my $MethodType = ($m or $method or "subroutine");
my $MinLength  = ($n or $nline or 10);
my $Exceptions = ($e or $except or $defaultexception);

# Allow in-place editing
$^I = "";

use strict;

my $SCRIPT  = "AddTest.pl";
my $ERROR   = "ERROR in $SCRIPT:";
my $WARNING = "WARNING in $SCRIPT:";

my @source = @ARGV;

if($Help or not @source){
    print "
This script can be used to add calls of test_start and test_stop
subroutines from BATL_test and improve the formatting of the code.

Usage:

      $SCRIPT [-h] [-u|-r] [-t=LEVEL] [-m=TYPE] FILE1 [FILE2 ...]

-h -help     Print this help message.

-e=PATTERN   Do not add tests to methods matching pattern.
-except=PATTERN Default is $defaultexception

-l=LEVEL     Insert call test_start into methods up to a contained
-level=LEVEL level of LEVEL. Top level is 0. Maximum is 2 in Fortran90.
             Default value is 1.

-m=TYPE      Add tests for methods described by TYPE. Possible values
-method=TYPE are 'subroutine', 'function', 'all' or 'none' 
             Only the first character matters. Default is 'subroutine'.

-n=NLINE     Add tests only to methods that are at least NLINE long.
-nline=NLINE Default is 10

-u -undo     Undo changes by copying the *_orig_ file back.

-r -redo     Redo processing by copying the _orig_ file back and
             then redo the processing.

FILE1 FILE2 ... 
      Source files to be checked. Note that the script should be run in 
      the directory where the source files are.

Examples:

Add tests at default level and fix code in all Fortran 90 files:

     ./$SCRIPT *.f90

Undo changes

    ./$SCRIPT -u *.f90

Redo processing one file and add tests to all functions and subroutines:

     ./$SCRIPT -r -level=2 -method=all -nline=1 some_code.f90 

Do not add tests, just fix code

     ./$SCRIPT -method=none *.f90

";
    exit;
}

# Check arguments
die "$ERROR invalid method type=$MethodType\n" if $MethodType !~ /^[asfn]/;
$MethodType =~ s/^a.*/all/;
$MethodType =~ s/^s.*/subroutine/;
$MethodType =~ s/^f.*/function/;
$MethodType =~ s/^n.*/none/;

# New test variable names
my @testvar = 
    ("test_start", "test_stop", "lVerbose", "StringTest", 
     "iTest", "jTest", "kTest", "iBlockTest", "iProcTest", "iVarTest", 
     "iDimTest", "xTest", "yTest", "zTest"), 

# Simple Fortran types with possible (len=..) and (kind=..) attributes:
my $SimpleType = '(real|integer|logical|character)(\s*\([^\)]+\))?\b';

# Obsolete Fortran types with *NUMBER, e.g. real*8 character*10
my $ObsoleteType = '(real|integer|logical|character)\s*\*\s*\d+';

# Derived Fortran type
my $DerivedType = 'type\s*\(\s*\w+\s*\)';

# Any Fortran Type
my $AnyType = "($SimpleType|$ObsoleteType|$DerivedType)";

my $source;
if($Undo or $Redo){
    foreach $source (@source){
	if(-f $source."_orig_"){
	    rename $source."_orig_", $source;
	}else{
	    print "$WARNING missing file $source"."_orig_\n";
	}
    }
}

# Done for undo
exit 0 if $Undo;

foreach $source (@source){

    next if $source =~ /$Exceptions/;

    # read source file into an array of lines
    open(FILE, $source);
    my @lines = <FILE>;
    close(FILE);

    my $iLevel = -1;   # number of program units started - 1
    my $indent;        # indentation at the beginning of the "unit"
    my $nindent;       # number of spaces in the indentation
    my $useteststart;  # true if "use ... test_start" was found.
    my $newtest;       # true if new test should be added
    my $oldtest;       # true if old test is already present
    my $addnamesub;    # true if namesub should be declared
    my $namesubline;   # declaration of NameSub = ...
    my $dotestline;    # declaration of DoTest logical
    my $separatorline; # !---- line
    my $addcalltest=1; # true if "call test_start" should be added
    my $unittype;      # program, subroutine or function
    my $unitname;      # name of program unit
    my $declaration;   # true inside declaration part
    my $implicitnone;  # true if "implicit none" was found at level 0
    my $usemodmain;    # true inside use ModMain... declaration
    my $testarguments; # arguments for test_start and test_stop
    my $removeoktest;  # true while removing old oktest code
    my $interface;     # true inside interface .. end interface
    my $istart0;       # first line of top level unit
    my $istart;        # first line of current method
    my $iseparator;    # index of !---- separator line (if present)
    my $separatorcheck;# check if original separator is correct

    my $i = -1;
    foreach (@lines){
	# line index
	$i++;

	# remove trailing spaces
	s/\s+\n/\n/;

	# Fix copyright message
	s/Michigan, portions used/Michigan,\n!  portions used/;
	$_ = "" if /This code is a copyright/ and /2002/;

	# Replace old variable names with new names
	s/\btest_string\b/StringTest/gi;
	s/\bItest\b/iTest/gi;
	s/\bJtest\b/jTest/gi;
	s/\bKtest\b/kTest/gi;
	s/\bXtest\b/xTest/gi;
	s/\bYtest\b/yTest/gi;
	s/\bZtest\b/zTest/gi;
	s/\bBLKtest\b/iBlockTest/gi;
	s/\bPROCtest\b/iProcTest/gi;
	s/\bVARtest\b/iVarTest/gi;
	s/\bDIMtest\b/iDimTest/gi;
	s/\biBLK\b/iBlock/gi;
	s/\bnBLK\b/MaxBlock/gi unless $source =~ /ModSize(_orig)?\.f90/;
	s/\boktest\b/DoTest/gi;
	s/\boktest_?me\b/DoTestMe/gi;

	# Fix capitalization errors
	s/^(\s*(end\s+)?)Module/$1module/i;

	# remove original !==== separator lines 
	$_ = '' if /^\s*\!\=+\!?$/;

	# fix comments !blabla --> ! blabla
	# except for the !INPUT ARGUMENTS: type Protex documentation
	s/(\s*\![\!\$]*)(\w)/$1 $2/ if not /^\s*\![A-Z\/ ]+:$/;

	# Ignore interface .. end interface
	$interface = 1 if /^\s*interface\b/i;
	if($interface){
	    $interface = 0 if /^\s*end interface\b/i;
	    next;
	}

	# Find start of subroutines and functions
	if(/^(\s*)(program|module|(recursive\s+)?subroutine|(recursive\s+)?($AnyType\s+)?function)\s+(\w+)/i){

	    my $actualindent = $1;
	    $unittype = lc($2);
	    $unitname = $+;

	    # Ignore "recursive" and type before functions and subroutines
	    $unittype =~ s/.*(function|subroutine)$/$1/;

	    # Increase level, set and check indentation
	    $iLevel++;
	    $nindent = 2*$iLevel;
	    $indent  =  " " x $nindent;

	    # remember first line index so we can check the length
	    $istart = $i;

	    # remember start line index of top level unit
	    # so we can declare test variables
	    $istart0 = $i if $iLevel == 0; 

	    print "$WARNING Incorrect indentation in $source at line $i: $_"
		if $actualindent ne $indent;

	    $declaration = 1;                  # start of declaration part
	    $implicitnone = 0 if $iLevel == 0; # implicit none is not yet found

	    # Decide if call test_start/stop should be added at this level
	    if($iLevel <= $TestLevel and $MethodType =~ /all|$unittype/ 
	       and $unittype !~ /module|program/
	       and $unitname !~ /$Exceptions/
	       and $source ne "ModBatsrusUtility.f90"){
		$newtest = 1;
	    }else{
		$newtest = 0;
	    }
	    # $oldtest will be set to 1 if code already contains testing
	    $oldtest = 0;

	    # Guess if this subroutine is called once per block or not
	    # If yes, call test_start/stop with iBlock argument.
	    if(/\biBlock\b/i){
		$testarguments = "NameSub, DoTest, iBlock";
	    }else{
		$testarguments = "NameSub, DoTest";
	    }
	    
	    # Construct lines for declaring DoTest, NameSub, and !----- separator
	    $dotestline    = "$indent  logical:: DoTest\n";
	    $namesubline   = "$indent  character(len=*), parameter:: ".
		"NameSub = '$unitname'\n";
	    $separatorline = "";
	    $separatorline = "$indent  !" . "-" x (76-$nindent) . "\n"
		unless $unittype eq "module";
	    $iseparator = 0;     # line index of separator
            $separatorcheck = 0; # check if original separator is correct

	    #print "indent=$indent, unittype=$unittype, unitname=$unitname"
	    #	.", newtest=$newtest\n";

	    next;
	}
	if($declaration){
	    if(/^\s+implicit\s+none/i){
		$implicitnone = 1 if $iLevel == 0;
		$_ = "" if $iLevel > 0;
		next;
	    }

	    # check if there is a !---------- separator line
	    $iseparator = $i if s/^\s+\!\-\-\-\-+\!?$//;
	    
	    # Skip empty lines and comments
	    next if /^$/ or s/^\s*\n/\n/ or /^\s*\!/;

	    if($separatorcheck++ > 1){
               # After iseparator is set, we should only get here once
               # More than once means that the !--- was improperly placed
               $iseparator = 0;
               $separatorcheck = 0;
            }

            # Fix ModSomething  ,only : --> ModSomething, ONLY:
            s/\bonly\s*:\b/ONLY:/i;
	    s/,ONLY/, ONLY/;
	    s/\s+(,\s+ONLY)/$1/;

            # remove iTest ... iVarTest from use ModMain
	    $usemodmain=1 if /^\s+use ModMain/i;
	    if($usemodmain){
		$usemodmain = 0 unless /\&$/;
		s/\blVerbose|(String|i|j|k|iBlock|iProc|iVar|iDim|x|y|z)Test\b\s*,//g;
		s/(,\s*)?(lVerbose|String|i|j|k|iBlock|iProc|iVar|iDim|x|y|z)Test\b//;
		s/,\s*$/\n/;

		# remove line if no variables are left in it
		next if s/^\s+use ModMain\s*,\s+ONLY:\s*$//i or s/^\s*\&$//;
                # remove continuation from previous line if this line is empty
		$lines[$i-1] =~ s/,\s*\&$// if s/^\s*$//;
		next;
	    }

	    $useteststart = 1 if /\btest_start\b/ and $iLevel == 0;

	    # Remove original NameSub declarations for sake of uniformity
	    if(/^\s+character.*parameter\s*::\s*NameSub\s*=/){
		$addnamesub = 1;
		$_ = "";
		next;
	    }
	    # Remove multiline NameSub
	    if(/^\s+NameSub\s*=/){
		$addnamesub = 1;
		$lines[$i-1] = "";
		$_ = "";
		next;
	    }

            # skip continuation lines
	    next if $lines[$i-1] =~ /\&\s*(\!.*)?$/;

	    # skip use statements
	    next if s/(\s*)use\b/$1use/i;

	    # Remove old declarations of DoTest and DoTestMe
	    if(/^\s+logical\s*::\s*DoTest(Me)?\b/){
		s/DoTest(Me)?\b\s*(=\s*\.(true|false)\.\s*)?,?\s*//g;
		$_ = '' if /^\s+logical\s*::\s*$/;
		$oldtest = 1 unless $unittype eq "module";
		next;
	    }

	    # If DoTest is an argument, no need to add test
	    $newtest = 0 if /^\s*logical.*intent.*DoTest/;

	    # skip all other variable declarations
	    #next if /^\s*($AnyType)/i;

	    if(/^[^\!]*\:\:/){
		next;
	    }

	    # Where to insert separator line and related lines?
	    if($iseparator){
		# Remove the original separator line to get length exact
		$lines[$iseparator] = '';
	    }else{
		# Add separator line to the previous line
		$iseparator = $i - 1;
	    }

	    # End of declarations 
	    $declaration = 0;
	}
	# fix code after the declaration part

        $addcalltest = 0 if /^\s+call test_start\b/;

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
	if(not /\!.*\b(cycle|exit|return|go\s*to)\b/){
	    s/\bcycle\b/CYCLE/;
	    s/\bexit\b/EXIT/;
	    s/\breturn\b/RETURN/;
	    s/\bgo\s*to\b/GOTO/;
	}

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
	if(/^(\s*)(contains|end\s+(program|module|subroutine|function))\b/){

	    my $addtest = ($oldtest or ($newtest and $i >= $istart+$MinLength));

	    if($Verbose and $unitname =~ /$Verbose/){
		print "!!! unitname=$unitname unittype=$unittype line $i: $_
addtest=$addtest oldtest=$oldtest newtest=$newtest
istart=$istart MinLength=$MinLength iseparator=$iseparator
";
	    }


	    # Declare DoTest if needed
	    $lines[$iseparator] .= $dotestline if $addtest;

	    # Add NameSub declaration if needed
	    $lines[$iseparator] .= $namesubline if $addnamesub or $addtest;
	    $addnamesub = 0;

	    # insert !----- separator line for subroutines and functions
	    $lines[$iseparator] .= $separatorline unless /^\s*end module/;
	    $separatorline = "";

	    # add call test_start line
	    $lines[$iseparator] .= "$indent  call test_start($testarguments)\n"
		if $addtest and $addcalltest;

	    # add call test_stop line
	    $_ = "$indent  call test_stop($testarguments)\n$_"
		if $addtest and $addcalltest;

	    if(not /\bcontains\b/){
		# Reduce level at end of program unit
		$iLevel--; 
		$nindent -= 2;
		$indent = " " x $nindent;
	    }

	    # Add !========= separator line and declare test stuff
	    if($iLevel >= 0){
		$_ .= "$indent  !" . ("=" x (76-$nindent)) . "\n";
	    }else{
		$_ .= "!" . ("=" x 78) . "\n";
		# Add declaration of test variables and methods 
		# to the beginning of the main unit if not yet done
		if(not $useteststart){
		    my $usetest;
		    my $text = join('', @lines[$istart0..$i]);
		    my $testvar;
		    foreach $testvar (@testvar){
			$usetest .= ", $testvar" if $text =~ /\b$testvar\b/i;
		    }
		    if($usetest){
			$usetest =~ s/^, //;
			#$usetest =~ s/$testvar$/ \&\n      $testvar"/
			#	    if length($usetext)>70;
			$lines[$istart0] .= "\n  use BATL_lib, ONLY: &\n       $usetest\n"
		    }
		}
	    }
		
	    #print "level=$iLevel, nindent=$nindent after $_";

	    # Reset some values
	    $unitname = "";
	    $newtest = 0;
	    $oldtest = 0;
	    $testarguments = "";
	}
    }

    my $text = join('', @lines);
    # remove double/triple empty lines
    $text =~ s/\n\n\n+/\n\n/g; 

    # Save copy of original
    rename $source, $source."_orig_";

    # write new file
    open(FILE, ">$source");
    print FILE $text;
    close FILE;
}
