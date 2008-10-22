#!/usr/bin/perl -s

my $Help    = ($h or $help);
my $Debug   = ($d or $debug);
my $Verbose = ($v or $verbose);
use strict;
&print_help if $Help;

#Get information for generating plots
print STDOUT "Enter Alt Slice For plotting: ";
my $AltSlice = <STDIN>;
chop($AltSlice);

print STDOUT "Enter number of plots: ";
my $nPlot = <STDIN>;
chop($nPlot);

print STDOUT "Enter function: ";
my $func = <STDIN>;
chop($func);
print "Function to plot: $func\n" if $Debug;

print STDOUT "Use autorange(y or n): ";
my $UseAutorange = <STDIN>;
chop($UseAutorange);
print "Using Autorange: $UseAutorange\n" if $Debug;

my $fmin = 0;
my $fmax = 1;
if ($UseAutorange eq 'n'){
    print STDOUT "Enter min(function): ";
    $fmin = <STDIN>;
    chop($fmin);
    
    print STDOUT "Enter max(function): ";
    $fmax = <STDIN>;
    chop($fmax);
}
my $MovieFile = "plot_movie_iAlt${AltSlice}.dat";

# Setup initial plotting commands for IDL
my @InitialCommands=("filename='$MovieFile'\n", "!x.range=[-.75,.75]\n", "!y.range=[-.75,.75]\n", "plotmode='contbargrid'\n","transform='n'\n", "func='$func'\n","npict=1\n","autorange='$UseAutorange'\n","fmax=$fmax\n","fmin=$fmin\n","setdevice,'plot_1.ps'\n",".r getpict\n",".r plotfunc\n",".r PwAxis\n",".r closedevice\n");

open(IDL, "| idl");

# Generate first plot
print IDL @InitialCommands;

# Loop over remaining snapshots, generating command to create each plot and 
# then giving command list to IDL for plotting
my $iSnapshot;
for $iSnapshot (2..$nPlot){
    my @Commands=("npict=$iSnapshot\n","setdevice,'plot_$iSnapshot.ps'\n",".r getpict\n",".r plotfunc\n",".r PwAxis\n",".r closedevice\n");
    print IDL @Commands;
}

close IDL;
exit;
##############################################################################

sub print_help{
    print "
Purpose: Use IDL to generate plots from file created with PostProcessIdl.pl. 
         Each plot is saved in plot_[iSnapshot].ps. Function syntax is 
         compatable with that described in chapter 5 of BATS-R-US manual

Usage: PlotsPw.pl [-v] [-d] [-h]

-h -help    - print help message and exit
-v -verbose - print verbose information
-d -debug   - print debugging information
";
    exit;
}
