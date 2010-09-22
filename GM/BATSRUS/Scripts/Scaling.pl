#!/usr/bin/perl -s -i

my $CompileCode  = $compile;
my $CreateRundir = $rundir;
my $SetupRun     = $setup;
my $SubmitRun    = $submit;
my $Help         = ($h or $help);

use strict;

if($Help){
    print "
Edit Scripts/Scaling.pl and Scripts/Scaling.job* files as necessary.
To compile the code, create run directories, setup PARAM.in and job files, 
and submit all the runs:

Scripts/Scaling.pl -compile -rundir -setup -submit

You can do these steps separately too.
The machine name is determined automatically.
The script stops if it is not recognized.
";
    exit;
}

my $Machine = `hostname`;
$Machine =~ s/\d*\n//;
print "Machine=$Machine\n";
my $IsHera;
my $IsPfe;
$IsHera = 1  if $Machine eq "hera";
$IsPfe =  1  if $Machine eq "pfe";

die "Unknown machine=$Machine\n" unless $IsHera or $IsPfe;

# Number of cores to run on
my @nCore;
@nCore  = (8,64,512,4096);

# Where to put the run directories
# my $Dir = "/nobackup/gtoth1/SCALING";
my $Dir = "SCALING";

my $JobScript = "Scripts/Scaling.job".$Machine;

if($CompileCode){
    `Config.pl -e=HdCrash -u=Crash -g=4,4,4,700,1`;
    `make CRASH PIDL`;
}

my $nCore;
if($CreateRundir){
    `mkdir -p $Dir`;
    foreach $nCore (@nCore){
	`make rundir DEFAULT_EXE=CRASH.exe`;
	`gunzip -c dataCRASH/input/hyades2d_1.3ns.out.gz > run/hyades2d_1.3ns.out`;
	`cp Param/CRASH/PARAM.in.hydro_scaling run/PARAM.in`;
	`cp $JobScript run/job`;
	`mv run $Dir/run_n$nCore`;
    }
}

if($SetupRun){
    foreach $nCore (@nCore){
	my $res = $nCore**(1/3);
	my $nRootX   = $res*40;
	my $nRootYZ  = $res*4;
	my $rundir = "$Dir/run_n$nCore";

	# Set grid size in PARAM.in
	@ARGV = ("$rundir/PARAM.in");
	while(<>){
	    s/\d+(\s+nRootBlockX)/$nRootX$1/;
	    s/\d+(\s+nRootBlock[YZ])/$nRootYZ$1/;
	    print;
	}

	# Set number of nodes/cores in job
	@ARGV = ("$rundir/job");
	if($IsPfe){
	    my $nNode = int($nCore/8+0.99);
	    while(<>){
		s/(\#PBS -l select)=\d+/$1=$nNode/;
		print;
	    }
	}elsif($IsHera){
	    my $nNode = int($nCore/16+0.99);
	    while(<>){
                s/(\#MSUB -l nodes)=\d+/$1=$nNode/;
		s/(run_n|srun \-n)\d+/$1$nCore/;
                print;
            }
	}
    }
}

if($SubmitRun){
    foreach $nCore (@nCore){
	my $rundir = "$Dir/run_n$nCore";
	`cd $rundir; qsub job`                   if $IsPfe;
	`cd $rundir; msub job | tail -1 > jobid` if $IsHera;
    }
}
