#!/usr/bin/perl
#$Id$

use strict;
use warnings;
use Fcntl qw(:flock SEEK_END);

#the structure of the array FileMask: [mask,rm,preplot,send]  
my @FileMask;

#the structure of the array FileList: [FileName,rm,preplot,send]
my @FileList;

#constants
use constant _TRUE_ => 0;
use constant _FALSE_ => 1;

#variables used for the parallel execution of the script 
my $nTotalThreads=1;
my $ThisThread=0;
my @childs;

my $host;
my $dir=".";
my $wait=_FALSE_;

my @TEST=("-np=1", "-host=/Volumes/data/EUROPA-LOCAL-RUN--TEMP", "-preplot=\'*.dat\'", "-dir=.");

#read the argument line
foreach (@ARGV) {
  
  if (/^-np=(.*)/i) {$nTotalThreads=$1; next};
  if (/^-host=(.*)/i) {$host=$1; next};
  if (/^-dir=(.*)/i) {$dir=$1; next};
  if (/^-wait/i) {$wait=_TRUE_; next;}
  
  if (/^-send=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_FALSE_, 'preplot'=>_FALSE_,'send'=>_TRUE_});next};
  if (/^-send-rm=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_TRUE_, 'preplot'=>_FALSE_,'send' =>_TRUE_});next};
  
  if (/^-preplot=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_TRUE_, 'preplot'=>_TRUE_,'send'=>_TRUE_});next};
  
  if (/^-rm=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_TRUE_, 'preplot'=>_FALSE_,'send'=>_FALSE_});next};
  
  if (/^-help/i) {
    print "Example:\n";
    print "./ampsDataProcessing.pl -np=1 -host=tower-left.engin.umich.edu:/Volumes/data/EUROPA-LOCAL-RUN--TEMP -preplot=\'*.sdf\' -dir=./\n"; 
    print "The argument line:\n";
    print "-help             -> print the list of the arguments\n";
    print "-wait   ->  the script will wait for the new files being generated untill it is killed by user\n";
    print "-np [number]      -> the number of threads\n\n";

    print "-send=[mask for file search, separated by ',']    -> files send to the remote host without preprocessing, the files are NOT removed\n";
    print "-send-rm=[mask for file search, separated by ','] -> files send to the remote host without preprocessing, the files are will be removed\n\n";
    
    print "-preplot=[mask for file search, separated by ',']    -> preplot the data files and send; after the files sent ther are removed\n";

    print "-rm=[mask for file search, separated by ',']   -> remove the files\n";
    
    print "-host   -> the name of the remote host and directory where the files will be copied\n";
    
    print "-dir   -> the directory where the output files will be created\n";
    exit;
  }
}



print "Starting post processing script\n";
print "Total threads:  $nTotalThreads\n";
 
do {
  if (-d $dir) {
    #create the file list
    for (my $i=0;$i<=$#FileMask;$i++) {
      my @mask;
      my @files;
      
      @mask=split(',',$FileMask[$i]{'mask'});
      foreach (@mask) {
        $_=~s/\'//g;
        $dir=~s/\'//g;
        
        my $t=$dir."/".$_;   
        @files=glob $t;
        
        foreach (@files) {
          push(@FileList,{'file'=>$_,'rm'=>$FileMask[$i]{'rm'}, 'preplot'=>$FileMask[$i]{'preplot'}, 'send'=>$FileMask[$i]{'send'}});
        }      
      }     
    }

    #start slave processes	
    if (@FileList) {
      if ($nTotalThreads==1) {
        ProcessDataFiles($nTotalThreads-1);
      }
      else {
        for (my $count=0;$count<$nTotalThreads;$count++) {
          my $pid = fork();
        
          if ($pid) {
            # parent
            #print "pid is $pid, parent $$\n";
            push(@childs, $pid);
          } elsif ($pid == 0) {
             # child
             ProcessDataFiles($count);
             exit 0;
          } else {
             die "couldnt fork: $!\n";
          } 
        }
       
        foreach (@childs) {
          my $tmp = waitpid($_, 0);
        }
      }
      
    }
  }
  
  if ($wait==_TRUE_) {
    sleep(120);
  }
}
while ($wait == _TRUE_);
 
print "Done.\n";

#=============================== Process Data Files =============================
sub ProcessDataFiles {
  my $ThisThread=$_[0];
       
  for (my $i=0;$i<=$#FileList;$i++) {
    if ($i%$nTotalThreads==$ThisThread) {
      #process the file
      my $fname=$FileList[$i]{'file'};
      my $rm=$FileList[$i]{'rm'};
      my $send=$FileList[$i]{'send'};
      my $preplot=$FileList[$i]{'preplot'};
      
      #determine whether the file is complete and can be processes
      if (-e $fname) {
        open(MYFILE,"< $fname") || next;
        
        if (! flock(MYFILE, LOCK_EX)) {
          close(MYFILE);
          next;
        }
        
        flock(MYFILE, LOCK_UN);
        close(MYFILE);
      }
      
      #preplot the data file
      if ((-e $fname) && ($preplot == _TRUE_)) {
        my $t=$fname;
        $t=~s/.dat$/.plt/;

        
#         system("preplot $fname");
        print "preplot $fname\n";
        `preplot $fname`;
        
        
        
        if (-e $t) {
#           system("rm -f $fname");
          print "rm -f $fname\n";
          `rm -f $fname`;
          
          $fname=$t;
        }
                
        $rm=_TRUE_;
        $send=_TRUE_;
      }
      
      #send the data file
      if ((-e $fname) && ($send == _TRUE_)) {
#        system("scp $fname $host");
        print "scp $fname $host\n";
        `scp $fname $host`;
      }
      
      #remove the data file
      if ((-e $fname) && ($rm == _TRUE_)) {
#        system("rm -f $fname");
        print "rm -f $fname\n";
        `rm -f $fname`;
      }      
    }
  }
}
 


