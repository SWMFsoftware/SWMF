#!/usr/bin/perl -s
#
# Script downloads the magnetograms and processes them
#
# MDI.pl -r=nnnn (-t=mmmm) (-e) (-mpi='mpirun -np 2')
###########################################################
#  IDL program
#

my $CR = ($r);
my $Till = ($t or $CR);
my $DoExec = ($e);
my $MpiRun = ($mpi);
print " Carrington rotation: first CR= $CR\n";
print " Carrington rotation: last  $Till\n";

system("make HARMONICS") if $DoExec;

$IDL_PRO="fits_to_asciicr";
my $Executable="../../../bin/HARMONICS.exe";
print ("$MpiRun"."$Executable\n") if $MpiRun;
############################################################
# The carringon rotation cycle will start with this number +1
#
 
#$CR=1964;
#
#############################################################

#
$CR = $CR -1;
  until($CR==$Till)
  {
$CR=$CR+1;
#
###########################################################
# The url of the data, data location directory and filenames
#
# Start of the cycle on CR numbers
print " CR= $CR\n";
#
  $ftpsite = "soi.stanford.edu";
  $ftpdir  = "magnetic/synoptic/carrot/M/".$CR."/";
  $ftpname = "synop_Mr_0.".$CR.".fits";
#
############################################################
# Check if the magnetogram has been already downloaded    
#
      if (-e $ftpname)
      {
        print("The Magnetogram for CR=$CR already exists in the directory \n");   
          
      }  
#
      else 
      { 
#
############################################################
#  Getting the magnetogram, renaming,running the IDL program,
#  renaming the output files
#
         system("wget http://".$ftpsite."/".$ftpdir.$ftpname);
         system("cp *".$CR.".fits fitsfile.fits");
         system("echo $IDL_PRO, \"CR=\"$CR | idl");

	 system("$MpiRun"."$Executable") if $DoExec;

         system("mv fitsfile.dat fitsfile_".$CR.".dat ");
         system("mv fitsfile.H fitsfile_".$CR.".H ");
         system("mv fitsfile_tec.dat fitsfile_tec".$CR.".dat ");

       }
#################################################################
# The end of the cycle on CR number
#

  }
