#!/usr/bin/perl -s
my $Help    = ($h or $help);
my $Verbose = ($v or $verbose);
my $AbsTol  = ($a or $abs or 1e-30);
my $RelTol  = ($r or $rel or 1e-12);
my $MaxErr  = ($m or $max or 100);

use strict;

&print_help if $Help or not @ARGV;

my $WARNING="WARNING in DiffNum.pl";
my $ERROR="ERROR in DiffNum.pl";

die "$ERROR: there should be two file arguments!\n" unless $#ARGV == 1;

my $File1 = $ARGV[0];
my $File2 = $ARGV[1];

open(FILE1, $File1) or die "$ERROR: could not open $File1\n";
open(FILE2, $File2) or die "$ERROR: could not open $File2\n";

# Pattern for numbers

my $Number = '[\+\-]?\d+(\.\d*)?([deDE][\+\-]?\d+)?';

my $iLine1;
my $iLine2;
my $eof1;
my $eof2;
my $OrigLine1;
my $OrigLine2;
my $Line1;
my $Line2;
my $Found1;
my $Found2;
my $Num1;
my $Num2;
my $Error=0;

while($Error < $MaxErr){

  SEARCH: {
      $Found1 = 0;      
      if( $Line1 =~ s/($Number)//){
	  $Num1   = $1;
	  $Found1 = 1;
	  last SEARCH;
      }
      $Line1 = <FILE1> or last SEARCH;
      $OrigLine1 = $Line1;
      $iLine1++;
      redo SEARCH;
  }
    
  SEARCH: {
      $Found2 = 0;
      if( $Line2 =~ s/($Number)//){
	  $Num2   = $1;
	  $Found2 = 1;
	  last SEARCH;
      }
      $Line2 = <FILE2> or last SEARCH;
      $OrigLine2 = $Line2;
      $iLine2++;
      redo SEARCH;
  }

    print "Num1 = $Num1\n" if $Verbose and $Found1;
    print "Num2 = $Num2\n" if $Verbose and $Found2;

    last unless $Found1 and $Found2;

    # replace Fortran style D exponent with Perl's 'E'.
    $Num1 =~ s/[dD]/E/;;
    $Num2 =~ s/[dD]/E/;;

    my $Diff = abs($Num1 - $Num2);
    my $Avrg = 0.5*(abs($Num1)+abs($Num2));

    if( $Diff > $AbsTol and $Diff > $RelTol*$Avrg ){
	print "abs($Num1 - $Num2)=$Diff at lines $iLine1 and $iLine2:\n";
	print "< $OrigLine1";
	print "> $OrigLine2";
	$Error++;
    }
}

warn "There are extra numbers in $File1\n" if $Found2 and not $Found1;
warn "There are extra numbers in $File2\n" if $Found1 and not $Found2;

die "$ERROR: number of significant differences is $Error\n" if $Error;

exit 0;

##############################################################################
sub print_help{

    print "
Purpose:
    Compare two files and show significant numerical differences only. 
    The text between numbers are ignored. Differences in number formats
    are also ignored. Absolute and relative differences below a given
    threshold are also ignored.

Usage: DiffNum.pl [-a=ABSERR] [-r=RELERR] [-m=MAXERR] FILE1 FILE2

   -a=ABSERR     - ignore differences smaller than ABSERR.
                   Default value is 1e-30

   -r=RELERR     - ignore relative differences smaller than RELERR.
                   The relative difference is the 2*|a-b|/(|a|+|b|).
                   Default value is 1e-12

   -m=MAXERR     - Quit if the number of differences reaches MAXERR.
                   Default value is 100

Examples:

   Compare File1 and File2. Show differences with absolute value exceeding 
   1e-10 and relative value exceeding 1e-6. At most 10 differences are shown:

DiffNum.pl -a=1e-10 -r=1e-6 -m=10 File1 File2

";
    exit 0;

}

