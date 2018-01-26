#!/usr/bin/perl

use strict;
use warnings;

my ($Stat,$Date,$file);
my $WorkDir = "/nobackup/vtenishe/Tmp_AMPS_test"; 
chdir "$WorkDir"; 

print "Work Directory=$WorkDir\n";

opendir(DIR, $WorkDir) or die "Cannnot open $WorkDir";

my @files;
@files=readdir DIR ;
closedir DIR;

open (FILE,'>',"test.log") || die "Cannot open the log file";
print FILE "List of the jobs to be executed:\n";

foreach $file (@files) {
  next unless ($file =~ m/^test_amps.pleiades.all/);
  print FILE "$file\n";
}

close (FILE);

foreach $file (@files) {
   next unless ($file =~ m/^test_amps.pleiades.all/);

   #Submit the job file and wait untill the execution is complete
   open (FILE,'>>',"test.log") || die "Cannot open the log file";
   print FILE "Submitting job $file\n";
   close (FILE);    

   `rm -rf AmpsTestComplete`; 

   for (my $iTest=0;$iTest<2;$iTest++) {
    sleep 120; 
    `/PBS/bin/qsub $file 1>> test.log 2>\&1`; 

    sleep 30;
    $Stat=`qstat | grep vtenishe`;
    chomp $Stat;

    open (FILE,'>>',"test.log") || die "Cannot open the log file";
    print FILE "Job $file has been submitted @ ".`date`;
    print FILE "$Stat @ ".`date`;
 
    while ($Stat =~ m/AMPS_pfe/) {
      sleep 60;
      $Stat=`qstat | grep vtenishe`;
      chomp $Stat;

      print FILE "$Stat @ ".`date`;
    }

    close (FILE);

    #the execute of the job is completed
    #if "AmpsTestComplete" is generated -> the job completed sucessfully 

    if (-e 'AmpsTestComplete') {
      last;
    }
  }

   open (FILE,'>>',"test.log") || die "Cannot open the log file";
   print FILE "Job $file is completed @ \n".`date`;
   close (FILE);
}
    


 
