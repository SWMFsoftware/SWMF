#!/usr/bin/perl -s

my $Help    = ($h or $H or $help);
my $Verbose = ($v or $verbose);
my $Gzip    = ($g or $gzip);
my $Repeat  = ($r or $repeat);
my $Output  = ($o or $output);

use strict;

my $ERROR = "ERROR in PostProc.pl";
my $WARNING = "WARNING in PostProc.pl";

my $NameOutput;
if($Output){
    die "$ERROR: options -o(utput) and -r(epeat) cannot be combined!\n"
	if $Repeat;
    die "$ERROR: option -o(utput) requires the name of a directory!\n" 
	unless $#ARGV == 0;
    $NameOutput = $ARGV[0];
    die "$ERROR: directory or file $NameOutput already exists!\n"
	if -e $NameOutput;
    mkdir ($NameOutput, 0777)
	or die "$ERROR: could not mkdir $NameOutput\n";
}

&print_help if $Help;

my $Pwd = `pwd`; chop $Pwd;

# Name of the plot directories for various components
my %PlotDir = (
    "GM" => "GM/IO2",
    "IE" => "IE/ionosphere",
    "IH" => "IH/IO2",
    "IM" => "IM/plots",
    "SC" => "SC/IO2",
    "UA" => "UA/Output,UA/data"
	    );

REPEAT:{
    my $Dir;
    foreach $Dir (sort keys %PlotDir){
	next unless -d $Dir;

	my $PlotDir = $PlotDir{$Dir};

	# Find the actual plot directory
	if($PlotDir =~ /,/){
	    my @PlotDir;
	    @PlotDir = split(/,/,$PlotDir);
	    foreach (@PlotDir){
		if(-d $_){
		    $PlotDir{$Dir} = $_;
		    $PlotDir       = $_;
		    last;
		}
	    }
	}

	warn "$WARNING: plot directory $PlotDir is missing\n" 
	    unless -d $PlotDir;
	next unless -d $PlotDir;

	print "cd $Dir\n" if $Verbose;
	chdir $Dir 
	    or die "$ERROR: could not change directory to $Dir\n";

	# Post process files if necessary
	if($Dir eq "IE"){
	    if($Gzip){
		&shell("./pION -g");
	    }else{
		&shell("./pION");
	    }
	}elsif( $Dir =~ /^SC|IH|GM$/ ){
	    &shell("./pIDL");
	    if($Gzip){
		&shell("./pTEC g");
	    }else{
		&shell("./pTEC p r");
	    }
	}
	chdir $Pwd;
    }

    if($Repeat){
	sleep $Repeat;
	redo REPEAT;
    }
}

# Done except for collecting output files
exit 0 unless $Output;

# Collect plot directories into $NameOutput 
# and make empty plot directories if requested
my $Dir;
foreach $Dir (sort keys %PlotDir){
    next unless -d $Dir;
    my $PlotDir = $PlotDir{$Dir};
    next unless -d $PlotDir;

    # Check if the plot directory is empty
    my @Files;
    opendir(DIR, $PlotDir)
	or die "$ERROR: could not open directory $PlotDir!\n";
    @Files = readdir(DIR) 
	or die "$ERROR: could not read directory $PlotDir!\n";
    closedir(DIR);
    if($#Files > 1){
	print "PostProc.pl: mv $PlotDir $NameOutput/$Dir with ",
	       $#Files-1," file"; print "s" if $#Files > 2; print "\n";
	rename $PlotDir, "$NameOutput/$Dir" or 
	    die "$ERROR: could not rename $PlotDir $NameOutput/$Dir\n";
	mkdir $PlotDir, 0777 or
	    die "$ERROR: could not mkdir $PlotDir\n";
    }else{
	warn "$WARNING: no files were found in $PlotDir\n";
    }
}

exit 0;

#############################################################

sub shell{
    my $command = join(" ",@_);
    print "$command\n" if $Verbose;
    my $result = `$command`;
    print $result if $Verbose;
}

##############################################################################
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Post-process plot files with Scripts/PostProc.pl}
#!ROUTINE: PostProc.pl - post-process plot files of the components
#!DESCRIPTION:
# This script is copied into the run directory and it should be executed there.
# The script post processes the plot files created by the components.
#
#!REVISION HISTORY:
# 02/12/2005 G. Toth - initial version
# 05/08/2005           added -o option to collect output into a directory tree
#EOP

sub print_help{
    print
#BOC
'Purpose:

   Post-process the plot files.

Usage:

   PostProc.pl [-h] [-v] [-g] [-r=REPEAT | -o DIR]

   -h -help    Print help message and exit.

   -v -verbose Print verbose information.

   -g -gzip    Gzip the ASCII files.

   -r=REPEAT   Repeat post processing every REPEAT seconds.
               Cannot be used with the -o option.

   -o -output  Output should be collected into a directory tree
               Cannot be used with the -r option. Requires DIR argument.

   DIR         Name of the directory tree for the -o option.
               The directory should be new.

Examples:

   Post-process the plot files:

PostProc.pl

   Post-process the plot files, compress them and print verbose info:

PostProc.pl -g -v

   Repeat post-processing every 360 seconds, gzip files and pipe 
   standard output and error into a log file:

PostProc.pl -r=360 -g >& PostProc.log &

   Collect processed output into a directory tree named OUTPUT/New:

PostProc.pl -o OUTPUT/New'

#EOC
    ,"\n\n";
    exit;
}
##############################################################################
