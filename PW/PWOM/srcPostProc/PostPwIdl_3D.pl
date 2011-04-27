#!/usr/bin/perl -s

my $Help    = ($h or $help);
my $Binary  = ($b or $binary);
my $Debug   = ($d or $debug);
my $Verbose = ($v or $verbose);
my $Append  = ($M);
use strict;

&print_help if $Help;

my $plot;
my $file;
my %plot_data;
my @test;
my $i=1;
my $iLine;

my $Pwd     = `pwd`; chop $Pwd;
my $PostPwIdl = "$Pwd/PostPwIdl_3D.exe";

# Get parameters from STDIN
print STDOUT "Enter number of Time slices: ";
my $nTime = <STDIN>;
chop($nTime);

# Get line plot list
my @plot_file = sort glob "plots/north_plots_iline????.out";
my $nLine = @plot_file;

if(@plot_file){
    die "ERROR: $PostPwIdl is not available, please make PwIDL\n"
        unless -x $PostPwIdl;
}else{
    die "WARNING: no line files were found\n";
}

#open program
my $pid = open(PostProc,"| $PostPwIdl");

# Pass parameters
if ($Binary) {
    print PostProc "real8\n";
}else{
    print PostProc "ascii\n";
}
print PostProc "$nLine\n";
print PostProc "$nTime\n";

# Pass file names

for $file (@plot_file) {
    print PostProc "$file\n";
#    my $Output= `$PostPwIdl < $file 2>&1`;
}
close PostProc;

###############################
# Do Southern Hemisphere
###############################

# Get plot list and read data into hash array
my @plot_file = sort glob "plots/south_plots_iline????.out";
my $nLine = @plot_file;
if($nLine == 0){
    die("No Plots for Southern Hemisphere. Only slices generated for North.\n");
}

if(@plot_file){
    die "ERROR: $PostPwIdl is not available, please make PwIDL\n"
        unless -x $PostPwIdl;
}else{
    die "WARNING: no line files were found\n";
}

#open program
my $pid = open(PostProc,"| $PostPwIdl");

# Pass parameters
if ($Binary) {
    print PostProc "real8\n";
}else{
    print PostProc "ascii\n";
}
print PostProc "$nLine\n";
print PostProc "$nTime\n";

# Pass file names

for $file (@plot_file) {
    print PostProc "$file\n";
#    my $Output= `$PostPwIdl < $file 2>&1`;
}
close PostProc;


exit;
##############################################################################

sub print_help{
    print "
Purpose: Extract 2D slices in altitude for ploting with idl. This routine 
         reads all of the individual field-line output and extracts the 
         values at a given altitude. It is assumed that the individual lines 
         are in the plots/ directory

Usage: PostPwIdl.pl [-v] [-d] [-h] [-b]

-h -help    - print help message and exit
-b -binary  - tells code that the field lines are saved in binary
-v -verbose - print verbose information
-d -debug   - print debugging information

";
    exit;
}
