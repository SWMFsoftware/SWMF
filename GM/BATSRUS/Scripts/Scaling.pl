#!/usr/bin/perl -s -i

my $Help         = ($h or $help);
&print_help if $Help;

my $Verbose      = $v;
my $Dryrun       = $n;
my $WeakScaling  = $weak;
my $CompileCode  = $compile;
my $CreateRundir = $rundir;
my $SubmitRun    = $submit;

use strict;

my $Machine = `hostname`;
$Machine =~ s/\d*\n//;
print "Machine=$Machine\n";
my $IsHera;
my $IsPfe;
$IsHera = 1  if $Machine eq "hera";
$IsPfe =  1  if $Machine eq "pfe";

die "Unknown machine=$Machine\n" unless $IsHera or $IsPfe;

# Number of nodes and cores to run on
my $nNode;
my $nCore;
my @nCore;
if($WeakScaling){
    @nCore  = (8,64,512,4096);
}else{
    @nCore  = (128,256,512,1024,2048,4096,8192);
}

# Input file
my $HyadesFile = "hyades2d_1.3ns.out";

# Where to put the run directories
# my $Dir = "/nobackup/gtoth1/SCALING";
my $Dir = "SCALING";
&shell("mkdir -p $Dir");
my $Rundir = "SCALING/run";

my $JobScript = "Scripts/Scaling.job.".$Machine;

if(not $WeakScaling){

    # Strong scaling uses a single run directory with many 
    # executables, plot directories and job scripts
    if($CreateRundir){
	&make_rundir;
	# Modify PARAM.in
	if($Dryrun){
	    print "edit $Rundir/PARAM.in\n";
	}else{
	    @ARGV = ("$Rundir/PARAM.in");
	    while(<>){
		s/^PLOTDIR/\#PLOTDIR/;  # use separate plot directories
		s/^AMR/\#AMR/;          # switch on AMR
		print;
	    }
	}

	foreach $nCore (@nCore){
	    $nNode = $nCore/8  if $IsPfe;
	    $nNode = $nCore/16 if $IsHera;
	    &shell("cp $JobScript $Rundir/job_$nCore");
	    &edit_jobscript("$Rundir/job_$nCore");
	    &shell("mkdir $Rundir/plot_$nCore/");
	}
	print "Creating $Rundir done\n";
    }
    if($CompileCode){
	&shell("Config.pl -e=HdCrash -u=Crash -nMaterial=3");
	foreach $nCore (@nCore){
	    my $nBlock = int(2*196608/$nCore + 0.99);
	    print "Compiling $Rundir/CRASH_$nCore.exe for nBlock=$nBlock\n";
	    &shell("Config.pl -g=4,4,4,$nBlock,1");
	    &shell("make CRASH");
	    &shell("cp src/CRASH.exe $Rundir/CRASH_$nCore.exe");
	}
	print "Compilations done\n";
    }

    if($SubmitRun){
	foreach $nCore (@nCore){
	    &submit_run($Rundir,"job_$nCore");
        }
    }

}else{
    # Weak scaling uses many run directories and a single executable
    if($CompileCode){
	&shell("Config.pl -e=HdCrash -u=Crash -g=4,4,4,700,1");
	&shell("make CRASH");
    }

    if($CreateRundir){
	foreach $nCore (@nCore){
	    $Rundir = "$Dir/run_n$nCore";
	    &make_rundir;
	    &shell("cp $JobScript $Rundir/job");

	    my $res = $nCore**(1/3);
	    my $nRootX   = $res*40;
	    my $nRootYZ  = $res*4;
	    my $rundir = "$Dir/run_n$nCore";

	    if($Dryrun){
		print "edit $Rundir/PARAM.in\n";
	    }else{
		# Set grid size in PARAM.in
		@ARGV = ("$Rundir/PARAM.in");
		while(<>){
		    s/\d+(\s+nRootBlockX)/$nRootX$1/;
		    s/\d+(\s+nRootBlock[YZ])/$nRootYZ$1/;
		    print;
		}
	    }
	    # Set number of nodes/cores in job
	    &edit_jobscript("$rundir/job");
	}
    }

    if($SubmitRun){
	foreach $nCore (@nCore){
	    &submit_run("$Dir/run_n$nCore","job");
	}
    }
}

exit;
#################################################################################
sub make_rundir{
    die "$Rundir already exists\n" if -d $Rundir and not $Dryrun;
    &shell("make rundir RUNDIR=$Rundir DEFAULT_EXE=CRASH.exe");
    &shell("gunzip -c dataCRASH/input/$HyadesFile.gz > $Rundir/$HyadesFile");
    &shell("cp Param/CRASH/PARAM.in.hydro_scaling $Rundir/PARAM.in");
}
#################################################################################
sub submit_run{
    my $rundir = shift;
    my $job    = shift;
    &shell("cd $rundir; qsub $job") if $IsPfe;
    &shell("cd $rundir; msub $job | tail -1 > ${job}id") if $IsHera;
}
#################################################################################
sub edit_jobscript{
    my $jobscript = shift;

    if($Dryrun){
	print "edit $jobscript\n";
	return;
    }

    @ARGV = ($jobscript);
    while(<>){
	if($IsPfe){
	    $nNode = int($nCore/8+0.99);
	    # Change number of nodes
	    s/(^\#PBS -l select)=\d+/$1=$nNode/;
	    s/(^\#PBS -q) normal/$1 wide/ if $nCore >= 1024;
	    s/(^\#PBS -q) wide/$1 normal/ if $nCore <  1024;
	}elsif($IsHera){
	    my $nNode = int($nCore/16+0.99);
	    s/(\#MSUB -l nodes)=\d+/$1=$nNode/;
	    s/(run_n|srun \-n)\d+/$1$nCore/;
	}
	if(not $WeakScaling){
	    # Change plot directory
	    s/plot_\d+/plot_$nCore/;
	    # Change executable and runlog filenames
	    s/CRASH\w*\.exe > runlog\w*/CRASH_$nCore.exe > runlog_$nCore/;
	}
	print;
    }
}
#################################################################################
sub shell{
    my $command = shift;
    print "$command\n" if $Verbose;
    `$command` unless $Dryrun;
}
#################################################################################
sub print_help{
    print "
This script can setup weak and strong scaling runs on hera and Pleiades.
It maybe necessary to edit Scripts/Scaling.pl and Scripts/Scaling.job* 
files to customize the code. 

Usage:

Scripts/Scaling.pl [-v] [-n] [-weak] [-compile] [-rundir] [-setup] [-submit]

-v        Switch on verbose mode.
-n        Switch on dry mode (do not execute commands).
-weak     Do weak scaling. Default is strong scaling.
-compile  Compile the executables.
-rundir   Create the run directory/directories.
-submit   Submit the jobs.

You can do these steps separately too.
The machine name is determined automatically.
The script stops if it is not recognized.
";
    exit;
}
