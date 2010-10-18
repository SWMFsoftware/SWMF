#!/usr/bin/perl -s -i

my $Help         = ($h or $help);
&print_help if $Help;

my $Verbose      = $v;
my $Dryrun       = $d;
my $WeakScaling  = $weak;
my $CompileCode  = $compile;
my $CreateRundir = $rundir;
my $SubmitRun    = $submit;
my $RadHydro     = $radhydro;
my $nCores       = $n;

use strict;

my $Machine = `hostname`;
$Machine =~ s/\d*\n//;
print "Machine=$Machine\n" if $Verbose;
my $IsHera;
my $IsPfe;
$IsHera = 1  if $Machine eq "hera";
$IsPfe =  1  if $Machine eq "pfe";

die "Unknown machine=$Machine\n" unless $IsHera or $IsPfe;

# Number of nodes and cores to run on
my $nNode;
my $nCore;
my @nCore;
if($nCores){
    @nCore  = split(/,/,$nCores);
}elsif($WeakScaling){
    @nCore  = (8,64,512,4096);
}else{
    @nCore  = (128,256,512,1024,2048,4096,8192);
}
print "Number of cores=@nCore\n" if $Verbose;

# Input files
my $ParamFile;
if($RadHydro){
    $ParamFile = "Param/CRASH/PARAM.in.scaling_radhydro";
}else{
    $ParamFile = "Param/CRASH/PARAM.in.scaling_hydro";
}
my $JobScript  = "Scripts/Scaling.job.".$Machine;
my $HyadesFile = "hyades2d_1.3ns.out";

# Where to put the run directories
my $Dir = "SCALING";
&shell("mkdir $Dir") unless -d $Dir;
my $Rundir = "SCALING/run"; # default for strong scaling

if($CompileCode){
    if($RadHydro){
	&shell("Config.pl -e=Crash -nMaterial=5 -nWave=30 -u=Crash");
    }else{
	&shell("Config.pl -e=HdCrash -nMaterial=3 -u=Crash");
    }
}

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

	# Create job scripts and plot directories
	foreach $nCore (@nCore){
	    &shell("cp $JobScript $Rundir/job_$nCore");
	    &edit_jobscript("$Rundir/job_$nCore",$nCore);
	    &shell("mkdir $Rundir/plot_$nCore/");
	}
	print "Creating $Rundir done\n";
    }
    if($CompileCode){
	foreach $nCore (@nCore){
	    my $nBlock;
	    if($RadHydro){
		$nBlock = int(128000/$nCore + 0.99);
		&shell("Config.pl -g=4,4,4,$nBlock,$nBlock");
	    }else{
		$nBlock = int(2*196608/$nCore + 0.99);
		&shell("Config.pl -g=4,4,4,$nBlock,1");
	    }
	    print "Compiling $Rundir/CRASH_$nCore.exe for nBlock=$nBlock\n";
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
	if($RadHydro){
	    &shell("Config.pl -g=8,8,8,100,100");
	}else{
	    &shell("Config.pl -g=8,8,8,100,1");
	}
	&shell("make CRASH");
    }

    if($CreateRundir){
	foreach $nCore (@nCore){
	    $Rundir = "$Dir/run_n$nCore";
	    &make_rundir;
	    &shell("cp $JobScript $Rundir/job");
	    &shell("cp src/CRASH.exe $Rundir/");

	    my $res = int($nCore**(1/3)+0.5);
	    my $nRootX   = $res*20;
	    my $nRootYZ  = $res*2;
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
	    &edit_jobscript("$rundir/job",$nCore);
	}
    }

    if($SubmitRun){
	foreach $nCore (@nCore){
	    &submit_run("$Dir/run_n$nCore","job");
	}
    }
}

exit;
###############################################################################
sub make_rundir{
    die "$Rundir already exists\n" if -d $Rundir and not $Dryrun;
    &shell("make rundir RUNDIR=$Rundir");
    &shell("gunzip -c dataCRASH/input/$HyadesFile.gz > $Rundir/$HyadesFile");
    &shell("cp $ParamFile $Rundir/PARAM.in");
}
###############################################################################
sub submit_run{
    my $rundir = shift;
    my $job    = shift;
    &shell("cd $rundir; qsub $job") if $IsPfe;
    &shell("cd $rundir; msub $job | tail -1 > ${job}id") if $IsHera;
}
###############################################################################
sub edit_jobscript{
    my $jobscript = shift;
    my $nCore     = shift;

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
	    s/\#+ (MSUB -l qos=exempt)/\#$1/ if $nCore > 4096;
	    s/(run_n|srun \-n)\d+/$1$nCore/;
	}
	if(not $WeakScaling){
	    # Change plot directory
	    s/plot_\d+/plot_$nCore/;
	    # Change executable and runlog filenames
	    s/CRASH\w*\.exe/CRASH_$nCore.exe/;
	    s/(runlog[^\d\n]*)\d*/$1_$nCore/;
        }
	print;
    }
}
###############################################################################
sub shell{
    my $command = shift;
    print "$command\n" if $Verbose;
    `$command` unless $Dryrun;
}
###############################################################################
sub print_help{
    print "
This script can setup weak and strong scaling runs on hera and Pleiades.
It maybe necessary to edit Scripts/Scaling.pl and Scripts/Scaling.job* 
files to customize the code. 

Usage:

Scripts/Scaling.pl [-v] [-d] [-n=CORES] [-weak | -radhydro] 
                   [-rundir] [-compile] [-submit]

-v        Switch on verbose mode.
-d        Dry run mode (do not execute commands).
-n=CORES  Set number of cores as a comma separated list of numbers.
          Default depends on scaling type and machine.
-weak     Do weak scaling. Default is strong scaling.
-radhydro Do 3D radhydro problem (only strong scaling). Default is 3D hydro.
-rundir   Create the run directory/directories (step 1)
-compile  Compile the executable(s) (step 2)
-submit   Submit the jobs (step 3)

Examples:

Create rundirectory for strong scaling of radhydro problem:

  Scripts/Scaling.pl -radhydro -rundir

Compile executables for 128 and 256 cores:

  Scripts/Scaling.pl -n=128,256 -radhydro -compile

Submit job for 128 cores:

  Scripts/Scaling.pl -n=128 -radhydro -submit

Show what would be done for a full weak scaling of the 3D hydro problem:

  Scripts/Scaling.pl -d -v -weak -rundir -compile -submit

";
    exit;
}
