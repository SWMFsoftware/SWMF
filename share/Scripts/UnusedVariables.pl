#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help = $h or $help;
my $Remove = $r;
my $Debug = $D;

# Allow in-place editing
$^I = "";

use strict;

my @source = @ARGV;

if($Help or not @source){
    print '
This script can be used to remove unused variables from Fortran code.
It uses the output of the NAG Fortran compiler to find the unused variables.
The script will not work with other compilers.
The -w compilation flag will be commented out in the CFLAG definition in 
Makefile.conf so that the compiler warnings about unused variables are kept.
It is a good idea to set zero level for optimization for sake of speed.

Usage:

      UnusedVariables.pl [-h] [-r] FILE1 [FILE2 ...]

-h    This help message

-r    Remove variables. Recompile code to check for errors.

FILE1 FILE2 ... 
      Source files to be checked. Note that the script should be run in 
      the directory where the source files are.

Examples:

See the unused variables for all Fortran files:

     UnusedVariables.pl *.f90 *.f

Remove unused variables from a certain file:

     UnusedVariables.pl -r amr.f90
';
    exit;
}

# Search for Makfile.conf up to 4 directories up
# Check if it is the NAG compiler
# Comment out the -w flag, so warnings are not suppressed
my $MakefileConf;
my $found;
my $compiler;
my $nag;
my $level;
foreach $level (0..4){
    $MakefileConf = "../" x $level . "Makefile.conf";
    next unless -s $MakefileConf > 100;
    $found=1;
    @ARGV = ($MakefileConf);
    while(<>){
	if(/FORTRAN_COMPILER_NAME=(.*)/){
	    $compiler = $1;
	    $nag = 1 if $compiler =~ /^(nagfor|f95)$/;
	}
	warn "Commenting out the -w flag in $MakefileConf\n" 
	    if $nag and s/^(CFLAG =.+) \-w(.+)/$1$2 # -w/;
	print;
    }
    last;
}

die "UnusedVariables.pl could not find Makefile.conf\n" unless $found;
die "UnusedVariables.pl only works with the nagfor compiler, not with $compiler\n"
    unless $nag;

my $source;
foreach $source (@source){

    next if $source eq "MpiTemplate.f90";

    `touch $source`;
    
    my $object = $source;
    $object =~ s/\.\w+$/.o/;

    my @log = `make $object 2>&1 1>/dev/null`;

    # no error messages
    next unless @log;

    # read source file into an array of lines
    open(FILE, $source);
    my @lines = <FILE>;
    close(FILE);

    foreach (@log){
	my $var;
	my $msg;

	print if $Debug;

	die "error compiling original $source: $_" if /^\s*(Fatal )?Error/;

	next unless s/^Warning: $source, line (\d+): //;
	my $nline = $1;

    	if (/(\w+) (explicitly imported) into/){
	    $var = $1;
	    $msg = $2;
	}elsif(/(Unused) (symbol|local variable|parameter) (\w+)/i){
	    next if /NAMEMOD|NAMESUB|DOTEST/; # Don't remove these
	    $msg = $1;
	    $var = $3;
	}elsif(/Local variable (\w+) is (initialised but never used)/){
	    $var = $1;
	    $msg = $2;
	}else{
	    next;
	}

	$msg  = lc($msg);
	my $line = $lines[$nline-1];
	my $method;
	if(/Unused parameter/i){
	    print "$source at line $nline: $var is $msg\n";
	}else{
	    $line =~ /end (module|function|subroutine) (\w+)/ or
		die "line=$line did not match end module/function/subroutine\n";
	    $method = "$1 $2";
	    print "$source at line $nline: $var is $msg in method=$method\n";
	}
	my $i;
	my $ilast = $nline;
	for($i = $nline-1; $i>0; $i--){
	    $line = $lines[$i-1];
	    last if $line =~ /$method/i; # stop at beginning of method
	    next if $line =~ /^\s*\!/;   # do not remove commented out code
	    next unless $line =~ /\b$var\b/i; # check if variable is present

	    $ilast = $i;
	}

	$line = $lines[$ilast-1];

	print "nline=$nline ilast=$ilast\n";

	next if $line =~ /IMPLEMENTED/; # do not remove these in ModUser

	print "original line $ilast:$line";

        # remove variable, variable = initial, variable(dimensions)
	$line =~ s/\b$var\b(\s*=[^,\n]+|\([^\)]+\))?//i; 

	$line =~ s/,\s*,/,/;             # remove double comma
	$line =~ s/:\s*,/:/;             # remove comma following :
	$line =~ s/,\s*\n/\n/;           # remove trailing comma
	$line =~ s/^(\s*),\s*/$1/;       # remove leading comma
	$line =~ s/^\s*\&\s*(\!.*)?\n//; # remove line "& !comment"

	# delete empty lines, empty declarations and use statements
	if($line =~ /(^\s*|::\s*|only\s*:\s*)(\!.*)?$/i){
	    $line = "";
	    # remove continuation from previous line
	    $lines[$ilast-2] =~ s/(\s*,\s*)?\&\s*\n/\n/;
	}

	print "modified line $ilast:$line";
	print "\n" if not $line;

	$lines[$ilast-1] = $line if $Remove;
    }

    if($Remove){
	my $orig = $source."_orig_";
	`mv $source $orig` unless -f $orig;
	open(FILE, ">$source");
	print FILE @lines;
	close FILE;

	unlink($object);
	my $log = `make $object 2>&1 1>/dev/null`;
	die "error compiling modified $source\n" if $log =~ /^\s*error/i;

    }
}
