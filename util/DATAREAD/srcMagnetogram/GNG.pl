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
my $DoExec = ($e or $mpi);
my $MpiRun = ($mpi);

print " Carrington rotation: first CR= $CR\n";
print " Carrington rotation: last  $Till\n";
print " The converted file is harmonics.dat\n" if $DoExec;

system("make HARMONICS") if $DoExec;

$IDL_PRO="fits_to_asciicr";
my $Executable="../../../bin/HARMONICS.exe";
$Executable=$MpiRun." ".$Executable if $MpiRun;

my $ftpsite = "gong2.nso.edu";
my $ftpdir = "QR/mqs/";

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

if (-e "CR".$CR."_GNG.fits")
      {
        print("The Magnetogram for CR=$CR is already downloaded \n");   
          
      }  
#
      else 
      { 

############################################################
#  Getting the magnetogram, renaming,running the IDL program,
#  renaming the output files
#


system("wget -r -nd -A \"*".$CR."_000.fits.gz\" ftp://".$ftpsite."/".$ftpdir);
system("gzip -d *".$CR."_000.fits.gz ");
system("mv *".$CR."_000.fits  CR".$CR."_GNG.fits");
system("rm *".$CR."_000.fits.gz");
      }

#
      if (-e "CR".$CR."_GNG.DAT")
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
	 system("cp CR".$CR."_GNG.fits fitsfile.fits");


system("echo $IDL_PRO, \"CR=\"$CR | idl");
system("mv fitsfile.dat CR".$CR."_GNG.DAT");
system("mv fitsfile.H fitsfile_".$CR.".H ");
system("mv fitsfile_tec.dat fitsfile_tec".$CR.".dat ");
      }
if ($DoExec){
    system("cp CR".$CR."_GNG.DAT  fitsfile.dat");
    system("$MpiRun "."$Executable");
    system("mv harmonics.dat CR".$CR."_GNG.dat");
    system("rm fitsfile.dat")
            }

#################################################################
# The end of the cycle on CR number
#

  }
