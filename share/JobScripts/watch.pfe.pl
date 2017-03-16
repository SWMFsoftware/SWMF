#!/usr/bin/perl
use strict;

my $pattern = $ARGV[0];

if(not $pattern){
    print "
Usage: watch.pfe.pl PATTERN >& watch.log &

Watch jobs matching pattern. If any of them starts to run, delete the others.
The following jobs are available now:

",`qstat -u \$USER`,"

Use PATTERN to select a subset of these, e.g.

watch.pfe.pl SWMF >& watch.log &
";
    exit;
}

my @results;
my $running;
my $user = $ENV{USER};
LOOP:{
    @results = `qstat -u \$USER | grep $pattern`;
    my $ids;
    foreach (@results){
	#/^(\d+)[^:]+:\d\d ([A-Z]) +((\dd\+)?\d\d:\d\d)/;
	/^(\d+).*([A-Z]) +(\S+)/;
	print "id=$1 status=$2 wait=$3\n";
	$ids .= " $1";
	$running = $1 if $2 eq "R";
    }
    print "-------------------------\n";
    if($running){
	$ids =~ s/ $running//;
	print "qdel $ids\n";
	`qdel $ids`;
	last LOOP;
    }

    sleep 5;
    redo;
}

print `qstat -u \$USER`;

print "Finished watch, job $running is running\n";
