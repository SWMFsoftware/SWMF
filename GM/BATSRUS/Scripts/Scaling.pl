#!/usr/bin/perl -s -i

my $CreateRundir = $create;
my $SetupRun     = $setup;
my $SubmitRun    = $submit;
my $Help         = ($h or $help);

use strict;

if($Help){
    print "
Edit Scripts/Scaling.pl and Scripts/Scaling.job* files as necessary.
To create run directories, setup PARAM.in and job files, and submit runs:

Scripts/Scaling.pl -create -setup -submit

You can do these steps separately too.
";
    exit;
}


# Number of cores to run on
my @nCore;
@nCore  = (8,64,512,4096);

# Where to put the run directories
my $Dir = "/nobackup/gtoth1/SCALING";
`mkdir -p $Dir`;

my $nCore;
if($CreateRundir){
    foreach $nCore (@nCore){
	`make rundir DEFAULT_EXE=CRASH.exe`;
	`gunzip -c dataCRASH/input/hyades2d_1.3ns.out.gz > run/hyades2d_1.3ns.out`;
	`cp Param/CRASH/PARAM.in.hydro_scaling run/PARAM.in`;
	`cp Scripts/Scaling.job.pfe run/job`;
	`mv run $Dir/run_n$nCore`;
    }
}

if($SetupRun){
    foreach $nCore (@nCore){
	my $res = $nCore**(1/3);
	my $nRootX   = $res*40;
	my $nRootYZ  = $res*4;
	my $rundir = "$Dir/run_n$nCore";
	@ARGV = ("$rundir/PARAM.in");
	while(<>){
	    s/\d+(\s+nRootBlockX)/$nRootX$1/;
	    s/\d+(\s+nRootBlock[YZ])/$nRootYZ$1/;
	    print;
	}

	`cp Scripts/Scaling.job.pfe $rundir/job`;
	my $nNode = $nCore/8;
	@ARGV = ("$rundir/job");
	while(<>){
	    s/(\#PBS -l select)=\d+/$1=$nNode/;
	    print;
	}
    }
}

if($SubmitRun){
    foreach $nCore (@nCore){
	my $rundir = "$Dir/run_n$nCore";
	`cd $rundir; qsub job`;
    }
}
