#!/usr/bin/perl
#$Id$
#create a copy of the original file that contains smaller number of lines; used in AMPS' test for comparing of the output files with the reference ones 

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);
use IO::File;


#arguments:
#$ARGV[0] -> the name of the data file
#$ARGV[1] -> step; only those lines which number is proportional to 'step' will be saved in the output  

open (fIn,"<$ARGV[0]");
open (fOut,">$ARGV[0].Reduced.Step=$ARGV[1]"); 

my $cnt=0;
my $step=$ARGV[1];
my $line;

while ($line=<fIn>) {
  if ($cnt%$step == 0) {
    print fOut "$line";
  }

  $cnt++;
} 

close (fIn);
close (fOut);

   






