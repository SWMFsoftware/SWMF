#!/usr/bin/perl -s

my $Help    = ($h or $H or $help);
my $Verbose = ($v or $verbose);
my $Gzip    = ($g or $gzip);
my $Repeat  = ($r or $repeat);
my $MakeMovie = ($m or $movie or $M or $MOVIE);
my $KeepMovieOnly = ($M or $MOVIE);
my $Rsync   = ($rsync or $sync);

use strict;

my $rsync = 'rsync -avz';
my $exclude = " --exclude '*.idl' --exclude '*.tec' --exclude '*.dat'".
    " --exclude '*.[hHTS]'";

my $INFO  = "PostProc.pl";
my $ERROR = "ERROR in PostProc.pl";
my $WARNING = "WARNING in PostProc.pl";

my $ParamIn = "PARAM.in";
my $RunLog  = "runlog runlog_[0-9]*";

my $NameOutput;
if(@ARGV){
    die "$ERROR: option -r(epeat) cannot be combined with output directory!\n"
	if $Repeat;
    die "$ERROR: only one directory name can be given!\n" 
	unless $#ARGV == 0;
    $NameOutput = $ARGV[0];
    die "$ERROR: directory or file $NameOutput already exists!\n"
	if -e $NameOutput;
    `mkdir -p $NameOutput`;
    die "$ERROR: could not mkdir -p $NameOutput\n" if $?;
}

die "$ERROR: option -rsync requires a target directory: rsync=TARGETDIR\n"
    if $Rsync eq "1";

&print_help if $Help;

my $Pwd = `pwd`; chop $Pwd;

# Set the movie option for pIDL
my $Movie;
if($KeepMovieOnly){
    $Movie = '-M';
}elsif($MakeMovie){
    $Movie = '-m';
}

# Name of the plot directories for various components
my %PlotDir = (
    "GM"     => "GM/IO2",
    "IE"     => "IE/ionosphere",
    "IH"     => "IH/IO2",
    "IM"     => "IM/plots",
    "PW"     => "PW/plots",
    "RB"     => "RB/plots",
    "SC"     => "SC/IO2",
    "UA"     => "UA/Output,UA/data",
    "STDOUT" => "STDOUT",
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
	    &shell("./pIDL $Movie");
	    if($Gzip){
		&shell("./pTEC A g");
	    }else{
		&shell("./pTEC A p r");
	    }
	}elsif( $Dir =~ /^IM/ ){
	    my @files=glob("plots/*.dat");
	    if($Gzip){
		&shell("gzip",@files) if @files;
	    }else{
		&shell("./Preplot.pl",@files) if @files;
	    }
	}elsif( $Dir =~ /^PW/ ){
	    if($Gzip){
		my @files=glob("plots/*.out");
		&shell("gzip",@files) if @files;
	    }
	}elsif( $Dir =~ /^RB/ ){
	    if($Gzip){
		my @files=glob("plots/*.fls");
		&shell("gzip",@files) if @files;
	    }
	}
	chdir $Pwd;
    }

    if($Rsync and not $NameOutput){
	my $Dir;
	foreach $Dir (keys %PlotDir){
	    my $PlotDir = $PlotDir{$Dir};
	    next unless -d $PlotDir;
	    my $command = $rsync;
	    $command .= $exclude if $Dir =~ /GM|SC|IH/;
	    &shell("$command $PlotDir/ $Rsync/$Dir") if -d $PlotDir;
	}
	&shell("$rsync $ParamIn $Rsync/")      if -f $ParamIn;
	&shell("$rsync runlog $Rsync/")        if -f "runlog";
	&shell("$rsync runlog_[0-9]* $Rsync/") if glob("runlog_[0-9]*");
    }

    if($Repeat){

	sleep $Repeat;
	redo REPEAT;
    }
}

# Done except for collecting output files
exit 0 unless $NameOutput;

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
	print "$INFO: mv $PlotDir $NameOutput/$Dir with ",
	       $#Files-1," file"; print "s" if $#Files > 2; print "\n";
	rename $PlotDir, "$NameOutput/$Dir" or 
	    die "$ERROR: could not rename $PlotDir $NameOutput/$Dir\n";
	mkdir $PlotDir, 0777 or
	    die "$ERROR: could not mkdir $PlotDir\n";
    }else{
	warn "$WARNING: no files were found in $PlotDir\n";
    }
}

# Copy and move some input and output files if present
if(-f $ParamIn){
    print "$INFO: cp $ParamIn $NameOutput/\n";
    `cp $ParamIn $NameOutput/`;
}else{
    warn "$WARNING: no $ParamIn file was found\n";
}
if(-f "runlog"){
    print "$INFO: mv runlog $NameOutput/\n";
    `mv runlog $NameOutput`;
}elsif(glob("runlog_[0-9]*")){
    print "$INFO: mv runlog_[0-9]* $NameOutput/\n";
    `mv runlog_[0-9]* $NameOutput`;
}else{
    warn "$WARNING: no $RunLog file was found\n";
}

print "$INFO: Restart.pl -o $NameOutput/RESTART\n";
&shell("Restart.pl -o $NameOutput/RESTART");

if($Rsync){
    print "$INFO: rsync -avz $NameOutput $Rsync\n";
    &shell("rsync -avz $NameOutput/ $Rsync");
    print "$INFO: rsync is complete\n";
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
#!QUOTE: \subsection{Post-Process Plot Files with PostProc.pl}
#!ROUTINE: PostProc.pl - post-process plot files of the components
#!DESCRIPTION:
# This script is copied into the run directory and it should be executed there.
# The script post processes the plot files created by the components.
# The script can run in the background and post process periodically.
# It can also collect the plot files, the standard output, the run log and 
# the PARAM.in file into an 'output directory tree'.
#
#!REVISION HISTORY:
# 02/12/2005 G. Toth - initial version
# 05/08/2005           added -o option to collect output into a directory tree
# 09/08/2005           for -o option copy PARAM.in and move runlog inot tree.
#EOP

sub print_help{
    print
#BOC
'Purpose:

   Post-process the plot files and/or collect them into an output tree.
   The PARAM.in and runlog files (if present) are also copied/moved into 
   the output tree.

Usage:

   PostProc.pl [-h] [-v] [-g] [-m | -M] [-r=REPEAT | DIR]

   -h -help    Print help message and exit.

   -v -verbose Print verbose information.

   -g -gzip    Gzip the ASCII files.

   -m -movie   Create movies from series of IDL files and keep IDL files.

   -M -MOVIE   Create movies from series of IDL files and remove IDL files.

   -r=REPEAT   Repeat post processing every REPEAT seconds.
               Cannot be used with the -o option.

   -rsync=TARGET Copy processed plot files into an other directory 
               (possibly on another machine) using rsync. The TARGET
               is the name of the target directory (with host machine). 
               When -rsync is used without the -o option, the original
               plot directories are synchronized. When -rsync is used with the
               -o option, the output directory is synchronized.
               rsync must be installed on the local and target machines,
               and no password should be required to execute rsync.

   DIR         Name of the directory tree to collect the processed files in.
               Cannot be used with the -r option. The directory is created
               and it should be new (to avoid overwriting older results).
               By default the processed data is not collected.

Examples:

   Post-process the plot files:

PostProc.pl

   Post-process the plot files and create movies from IDL output:

PostProc.pl -M

   Post-process the plot files, compress them, rsync to another machine
   and print verbose info:

PostProc.pl -g -rsync=ME@OTHERMACHINE:My/Results -v

   Repeat post-processing every 360 seconds, gzip files and pipe 
   standard output and error into a log file:

PostProc.pl -r=360 -g >& PostProc.log &

   Collect processed output into a directory tree named OUTPUT/New
   and rsync it into the run/OUTPUT/New directory on another machine:

PostProc.pl -rsync=ME@OTHERMACHINE:run/OUTPUT/New OUTPUT/New'

#EOC
    ,"\n\n";
    exit;
}
##############################################################################
