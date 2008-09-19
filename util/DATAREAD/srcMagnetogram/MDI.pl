#!/usr/bin/perl -s
#
# Script downloads the magnetograms and processes them
# r means Carrington rotation t - Till e - execute
# MDI.pl -r=nnnn (-t=mmmm) (-e) (-mpi='mpirun -np 2')
###########################################################
#  IDL program
#

my $CR = ($r);
my $Till = ($t or $CR);
my $DoExec = ($e or $mpi);
my $MpiRun = ($mpi);

print " Carrington rotation: first CR= $CR\n";
print " Carrington rotation: last  $Till\n";
print " The converted file is harmonics.dat\n" if $DoExec;

system("make HARMONICS") if $DoExec;

$IDL_PRO="fits_to_asciicr";
my $Executable="../../../bin/HARMONICS.exe";
$Executable=$MpiRun." ".$Executable if $MpiRun;

$ftpsite = "soi.stanford.edu";

print ("$Executable\n") if $DoExec;
############################################################
# The carringon rotation cycle will start with this number +1
#

$CR = $CR -1; 

#
#############################################################

#

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

#
############################################################
# Check if the magnetogram has been already downloaded    
#
      if (-e "CR".$CR."_MDI.fits")
      {
        print("The Magnetogram for CR=$CR is already downloaded \n");   
          
      }  
#
      else 
      { 
#
############################################################
#  Getting the magnetogram, renaming
  $ftpdir  = "magnetic/synoptic/carrot/M/".$CR."/";
  $ftpname = "synop_Mr_0.".$CR.".fits";

	system("wget http://".$ftpsite."/".$ftpdir.$ftpname);
        system("mv *".$CR.".fits CR".$CR."_MDI.fits");
      }


#
      if (-e "CR".$CR."_MDI.DAT")
      {
        print("The Magnetogram for CR=$CR is already converted \n");   
          
      }  
#
      else 
      {  
#
############################################################
#  Running the IDL program,
#  renaming the output files
#
	 system("cp CR".$CR."_MDI.fits fitsfile.fits");

         system("echo $IDL_PRO, \"CR=\"$CR | idl");
	 system("mv fitsfile.dat CR".$CR."_MDI.DAT ");
         system("mv fitsfile.H fitsfile_".$CR.".H ");
         system("mv fitsfile_tec.dat fitsfile_tec".$CR.".dat ");
      }
if ($DoExec){
    system("cp CR".$CR."_MDI.DAT  fitsfile.dat");
    system($Executable);
    system("mv harmonics.dat CR".$CR."_MDI.dat");
    system("rm fitsfile.dat")
            }

#################################################################
# The end of the cycle on CR number
#

  }
