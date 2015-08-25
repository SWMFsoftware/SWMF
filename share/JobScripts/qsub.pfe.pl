#!/usr/bin/perl
use strict;

my $script = shift(@ARGV);
my $name   = shift(@ARGV);
my @machine = @ARGV;

if(not $script or $script =~ /\-+h/i or not $name){
    print "
Usage: qsub.pfe.pl SCRIPT NAME [MACHINE1] [MACHINE2] [MACHINE3] ...

Submit generic job script to multiple machine types.
Use a unique NAME argument to identify the jobse.
Only the first four characters of the NAME are used.
If no machine is specified, 4 jobs will be submitted for the 4 machine
types (Westmere, IvyBridge, SandyBridge, Haswell). Otherwise,
the job will be submitted for the listed machines.
Only the first three characters of the machine types are used.
Use watch.pl to make sure that when any of the jobs start to run, 
the others get deleted from the queuewith qdel. Note you can
add or delete jobs with matching NAME while watch.pl is running.

Example:

qsub.pfe.pl job.long Mars
watch.pfe.pl Mars >& watch.log &
";
    exit;
}

# Read original script into $text
print "qsub.pfe.pl reading $script\n";
my $text;
open(SCRIPT, $script) or die "Could not open $script\n";
$text = join("", <SCRIPT>);
close SCRIPT;

# Copy original script
my $machine;
my @script;
@machine = ('Ivy', 'San', 'Has', 'Wes') if not @machine;

foreach $machine (@machine){
    $machine =~ s/(...).*/$1/;

    my $fileout = "$script.$machine";
    print "creating $fileout\n";

    # Change name of the job to show machine name
    $text =~ s/^(#PBS -N).*/$1 $name$machine/m;

    # Comment out all active machine selections
    $text =~ s/^#(PBS -l.*model=.*)$/### $1/m;
    
    # Uncomment the line for model=$machine
    $text =~ s/^### (PBS -l.*model=$machine)$/#$1/im;

    # Change the name of the resubmit script
    $text =~ s/^(if.*qsub) .*$/$1 $fileout/m;

    open(SCRIPT, ">$fileout") or die "Could not open $fileout\n";
    print SCRIPT $text;
    close SCRIPT;
}


# submit jobs;
foreach $machine (@machine){
    print "qsub $script.$machine\n";
    `qsub $script.$machine`;
}
