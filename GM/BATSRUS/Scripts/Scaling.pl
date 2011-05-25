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
my $BlockSize    = $b;
my $nCores       = $n;

use strict;

my $Machine = `hostname`;
$Machine =~ s/\d*\n//;
print "Machine=$Machine\n" if $Verbose;
my $IsHera;
my $IsPfe;
my $IsUbgl;

$IsHera = 1  if $Machine eq "hera";
$IsUbgl = 1  if $Machine eq "ubgl";
$IsPfe  = 1  if $Machine eq "pfe";

die "Unknown machine=$Machine\n" unless $IsHera or $IsPfe or $IsUbgl;

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

# Grid size
my $nBlock;
my $nImplBlock;

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

    $BlockSize=4 unless $BlockSize;
    my $nRootX  = 320/$BlockSize;
    my $nRootYZ = 32/$BlockSize;

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
		s/\#(CHECKGRIDSIZE)/$1/;# switch off CHECKGRIDSIZE
		s/\d+(\s+nRootBlockX)/$nRootX$1/;
		s/\d+(\s+nRootBlock[YZ])/$nRootYZ$1/;
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
	    if($RadHydro){
		$nBlock     = int(8192000/$BlockSize**3/$nCore + 0.99);
		$nImplBlock = $nBlock;
	    }else{
		$nBlock     = int(25165824/$BlockSize**3/$nCore + 0.99);
		$nImplBlock = 1;
	    }
	    &shell("Config.pl -g=$BlockSize,$BlockSize,$BlockSize,$nBlock,$nImplBlock");
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
    if($RadHydro){
	$BlockSize  = 4 unless $BlockSize;
	$nBlock     = 16384 / $BlockSize**3;
	$nImplBlock = $nBlock;
    }else{
	$BlockSize  = 8 unless $BlockSize;
	$nBlock     = 40960 / $BlockSize**3;
	$nImplBlock = 1;
    }

    if($CompileCode){
	&shell("Config.pl -g=$BlockSize,$BlockSize,$BlockSize,$nBlock,$nImplBlock");
	&shell("make CRASH");
    }

    if($CreateRundir){
	foreach $nCore (@nCore){
	    $Rundir = "$Dir/run_n$nCore";
	    &make_rundir;
	    &shell("cp $JobScript $Rundir/job");
	    &shell("cp src/CRASH.exe $Rundir/");

	    my $res = int($nCore**(1/3)*16/$BlockSize + 0.5);
	    my $nRootX   = $res*10;
	    my $nRootYZ  = $res;
	    my $rundir = "$Dir/run_n$nCore";

	    if($Dryrun){
		print "edit $Rundir/PARAM.in\n";
	    }else{
		# Set grid size in PARAM.in
		@ARGV = ("$Rundir/PARAM.in");
		while(<>){
		    s/\#(CHECKGRIDSIZE)/$1/;
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
    &shell("cd $rundir; msub $job | tail -1 > ${job}id") if $IsHera or $IsUbgl;
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
	}elsif($IsUbgl){
	    if($nCore <= 128 or $nCore == 512){
		# the minimum number of nodes is 128
		my $nNode = $nCore;
		$nNode = 128 if $nNode < 128;
		s/^(\#MSUB -l nodes)=\w+/$1=$nNode/;
		# run with 1 core per node using mpirun -np $nCore
		s/^mpirun (-np \d+|-mode VN)/mpirun -np $nCore/;
	    }else{
		# 2 cores / node and convert to kilonode notation
		my $nNode = int($nCore/2);
		$nNode = ($nNode/1024)."k" if $nNode >= 1024;
		s/^(\#MSUB -l nodes)=\w+/$1=$nNode/;
		# run with 2 cores per node, using mpirun -mode VN
		s/^mpirun (-np \d+|-mode VN)/mpirun -mode VN/;
	    }
	    s/^(\#MSUB -q) \w*/$1 pdebug/ if $nCore <= 1024;
	    s/^(\#MSUB -q) \w*/$1 pshort/ if $nCore >  1024 and $nCore <= 8192;
	    s/^(\#MSUB -q) \w*/$1 pbatch/ if $nCore >  8192;
	    my $rundir = `pwd`; chop $rundir; $rundir .= "/SCALING/run_n$nCore";
	    s/^cd .*/cd $rundir/;
	}
	if(not $WeakScaling){
	    # Change plot directory
	    s/\bplot\b/plot_$nCore/;
	    # Change executable and runlog filenames
	    s/CRASH\w*\.exe/CRASH_$nCore.exe/;
	    s/\b(runlog[^\d\n]*)\b/$1_$nCore/;
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
-b=BLOCKSIZE Set grid blocks to have BLOCKSIZE*BLOCKSIZE*BLOCKSIZE cells.
-n=CORES  Set number of cores as a comma separated list of numbers.
          Default depends on scaling type and machine.
-weak     Do weak scaling. Default is strong scaling.
-radhydro Do 3D radhydro problem. Default is 3D hydro.
-rundir   Create the run directory/directories (step 1)
-compile  Compile the executable(s) (step 2)
-submit   Submit the jobs (step 3)

Examples:

Compile executables for strong scaling radhydro on 128 and 256 cores:

  Scripts/Scaling.pl -n=128,256 -radhydro -compile

Create rundirectory for strong scaling of radhydro problem:

  Scripts/Scaling.pl -n=128,256 -radhydro -rundir

Submit job for 128 cores:

  Scripts/Scaling.pl -n=128 -radhydro -submit

Full weak scaling of the pure hydro problem with 16 cubed blocks:

  Scripts/Scaling.pl -b=16 -n=1,8,64,512,4096,32768 -weak -compile -rundir -submit

";
    exit;
}
