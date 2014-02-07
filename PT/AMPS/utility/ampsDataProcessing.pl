#!/usr/bin/perl
#$Id$

use strict;
use warnings;

my @ProcessDataFileMask=('GroundBasedObservation.*','ColumnDensityMap.spherical.*','NA.s=0','Sphere=0.NA.s=0','LimbColumnDensity*');
my @RemoveDataFileMask=('NAPLUS.s=1','Sphere=0.NAPLUS.s=1','distribution.NAPLUS.s=1','flux.NA.s=0','flux.NAPLUS.s=1');
my @CopyDataFileMask=('SurfaceProperties.*','distribution.NA.s=0','TailColumnDensity.*','LimbColumnDensity*','Moon.Anti-sunwardColumnIntegrals*','pic.Moon.Kaguya.TVIS*');
my $RemoteLocation="tower-right.engin.umich.edu:HEC"; 

use constant _REMOVE_DATA_FILE__ON_  => 0;
use constant _REMOVE_DATA_FILE__OFF_ => 1;

use constant _REMOTE_LOCATION_SEND_FILE__ON_  => 0;
use constant _REMOTE_LOCATION_SEND_FILE__OFF_ => 1;

use constant _WAIT_NEW_DATA_FILES__ON_  => 0;
use constant _WAIT_NEW_DATA_FILES__OFF_ => 1;


my $RemoveDataFileFlag=_REMOVE_DATA_FILE__OFF_;
my $RemoteLocationSendFlag=_REMOTE_LOCATION_SEND_FILE__OFF_;
my $WaitNewDataFileFlag=_WAIT_NEW_DATA_FILES__OFF_;


my $nTotalThreads=1;
my $ThisThread=0;
my @childs;

my @RemoveDataFileList;
my @CopyUnprocessedDataFileList;
my @PreplottedDataFileList;


#parse the argument line
#argument list
#-np [number]      -> the number of threads
#-send=[on,off]    -> send the files to the remote location
#-rm=[on,off]      -> remove post processed and send files

print "Starting post processing script\n";
print "Total threads:  $nTotalThreads\n";
 
do {
  #create the file lists
  #the list of files to be removed
  @RemoveDataFileList=();
	
  if ($RemoveDataFileFlag==_REMOVE_DATA_FILE__ON_) {
    foreach (@RemoveDataFileMask) {
	  push(@RemoveDataFileList,<pic.$_.out=*.dat*>);
	}
  }
	
  #copy unprocessed data files
  @CopyUnprocessedDataFileList=();
	 
  if ($RemoteLocationSendFlag==_REMOTE_LOCATION_SEND_FILE__ON_) {
	foreach (@CopyDataFileMask) {
      push(@CopyUnprocessedDataFileList,<pic.$_.out=*.dat*>);
    }
  }
  
  #the list of files to be pre-ploted 
  @PreplottedDataFileList=();
  
  foreach (@ProcessDataFileMask) {
    push(@PreplottedDataFileList,<pic.$_.out=*.dat*>);
  } 
	
  #start slave processes	
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
  
   if ($WaitNewDataFileFlag == _WAIT_NEW_DATA_FILES__ON_) {
   	 sleep(120);
   }
}
while ($WaitNewDataFileFlag == _WAIT_NEW_DATA_FILES__ON_);
 
print "Done.\n";

#=============================== Process Data Files =============================
sub ProcessDataFiles {
  my $ThisThread=shift;
  my $FileCounter=0;
  
  #remove unwanted files
  foreach (@RemoveDataFileList) {
    my $fname="$_";
    my $fn;

    if ($FileCounter%$nTotalThreads==$ThisThread) {
      if (open $fn,"+<",$fname) {
        close $fn;
               
        if ($RemoveDataFileFlag == _REMOVE_DATA_FILE__ON_) {       
          print "Remove $fname\n";
	      $fname=~s/\(/\\\(/g;
	      $fname=~s/\)/\\\)/g;
	      `rm -f $fname`;
        }
      }
    }
    
    $FileCounter++;
  } 
  
  #copy data files
  foreach (@CopyUnprocessedDataFileList) {
    my $fname="$_";
    my $fn;

    if ($FileCounter%$nTotalThreads==$ThisThread) {
      if (open $fn,"+<",$fname) {
        close $fn;
        
        $fname=~s/\(/\\\(/g;
        $fname=~s/\)/\\\)/g;
        
        print "Copy ./$fname to $RemoteLocation\n";
        
        if ($RemoteLocationSendFlag==_REMOTE_LOCATION_SEND_FILE__ON_) {
        	`scp ./$fname $RemoteLocation`;
        }
        
        
        if ($RemoveDataFileFlag == _REMOVE_DATA_FILE__ON_) {
          print "Remove $fname\n";
          `rm -f $fname`;
        }
      }
    }
    
    $FileCounter++;
  } 


  #preplot and scp other data files
  foreach (@PreplottedDataFileList) {
    my $fh;
    my $SourceFileName=$_;
    my $ProcessedFileName;
    my $FoundProcessedFile;

    if ($FileCounter%$nTotalThreads==$ThisThread) {
      if (open $fh,"+<",$SourceFileName) {
        close $fh;

        print "Thread $ThisThread, Preplot: $SourceFileName\n";
        `preplot $SourceFileName`;

         $FoundProcessedFile=0;

         #check if the processed file exists
         $ProcessedFileName=$SourceFileName;
         $ProcessedFileName=~s/\.dat$/\.\.plt/;

         # check if *..plt file is found -> rename the file         
         if (-e $ProcessedFileName) {           
           $FoundProcessedFile=1;
           
           my $newname=$ProcessedFileName; 
           $newname=~s/\.\./\./;               
           print "Thread $ThisThread, Rename $ProcessedFileName to $newname\n";    
           `mv $ProcessedFileName $newname`;      
           $ProcessedFileName=$newname;
         }
         
         if ($FoundProcessedFile==0) {
         	#check if *.plt file exists
         	
         	$ProcessedFileName=$SourceFileName;
         	$ProcessedFileName=~s/\.dat$/\.plt/;
         	
         	if (-e $ProcessedFileName) {
         	  $FoundProcessedFile=1;
         	}	
         }
         
         if ($FoundProcessedFile==0) {
           #if no post processed files are found => zip and copy the original data file

           `gzip $SourceFileName`;
           $ProcessedFileName=$SourceFileName.".gz";
           $FoundProcessedFile=1;

           print "Thread $ThisThread, Error: no processed file found for the source $SourceFileName, the name of *.zip file is $ProcessedFileName \n";
         }
         
         #copy the processed file
         if ($RemoteLocationSendFlag==_REMOTE_LOCATION_SEND_FILE__ON_) {
            print "Thread $ThisThread, Copy ./$ProcessedFileName to $RemoteLocation\n";
            `scp ./$ProcessedFileName $RemoteLocation`;
         }

         #remove unwanted files
         if ($FoundProcessedFile==1) {
           print "Thread $ThisThread, Remove $SourceFileName\n";
           $SourceFileName=~s/\(/\\\(/g;
           $SourceFileName=~s/\)/\\\)/g;
           `rm -f $SourceFileName`;         	
         }

         if ($RemoveDataFileFlag == _REMOVE_DATA_FILE__ON_) {     
           print "Thread $ThisThread, Remove $ProcessedFileName\n";
           $ProcessedFileName=~s/\(/\\\(/g;
           $ProcessedFileName=~s/\)/\\\)/g;
           `rm -f $ProcessedFileName`;
         }

       }
    }  
    
    $FileCounter++;
  }
}



