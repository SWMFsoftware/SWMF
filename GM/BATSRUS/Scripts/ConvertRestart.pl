#!/usr/bin/perl 
#^CFG COPYRIGHT UM

#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Change Endianness of BATSRUS Restart Files with Scripts/ConvertRestart.pl}
#!ROUTINE: ConvertRestart.pl - change the endianness of BATSRUS restart files
#!DESCRIPTION:
# When moving restart files from one machine to another, it may be
# necessary to change the endianness of the binary files.
# This script is specifically written to convert all the binary restart
# files produced by BATSRUS, as well as copy the ASCII header file.
#
#!REVISION HISTORY:
# 07/03/2001 G. Toth - initial revision
# 03/11/2003 G. Toth - generalized to work on a machine with 
#                      arbitrary endianness
#EOP
if($#ARGV != 1){
    print 
#BOC
"Purpose: change the endianness of blk*.rst and octree.rst files, 
         and copy restart.H.

This perl script can be run on a machine with arbitrary endianness,
however do not run this on the Cray, because the Perl interpreter 
does not interpret long integers correctly on the Cray. 

The typical usage is

    mkdir restart_export
    ConvertRestart.pl restartOUT restart_export
    tar cf restart_export.tar restart_export
    scp restart_export.tar OTHERMACHINE:BATSRUS/run

The files are read from the 'inputdir' (restartOUT) and the restult is 
saved into the 'outputdir' (restart_export).

If the files are moved to the Cray or from the Cray, the octree file needs
to be fixed for the 4 vs. 8 byte integers using FixI4toI8.pl or FixI8toI4.pl.
The endianness can be changed with ConvertRestart both after or before the 
8-byte-integer-fix."
#EOC
    ,"\n\n";
    exit;
}

$indir=$ARGV[0];
$outdir=$ARGV[1];

# Check things
opendir(INDIR,  $indir)  or die "Cannot open input directory $indir!\n";
opendir(OUTDIR, $outdir) or die "Cannot open output directory $outdir!\n";
closedir OUTDIR;

die "Input and output directories must be different!\n" if $indir eq $outdir;

# There is 1 octree.rst file and many blk files
@rstfiles = grep /^(blk\d*|octree)\.rst$/, readdir INDIR;
print "Number of blocks in $indir is ",$#rstfiles,"\n";

# No end of line in files
undef $/;

# Default integer format is VAX (little-endian, e.g. Linux PC)
$intform='V';

print "Converting indir=$indir to outdir=$outdir...\n";
foreach $rstfile (@rstfiles){

    open(INFILE, "$indir/$rstfile") or die "Cannot open $rstfile!\n";
    $_=<INFILE>;
    close INFILE;
    
  ENDIAN:{
      $len = unpack $intform,substr($_,0,4);

      if($len == 8 or $len == 12){
	  &convert4;
      }elsif($len == 16 or $len == 24){
	  &convert8;
      }elsif($intform eq 'V'){
	  # Try Network integer format (big-endian, e.g. SGI, IBM, Cray)
	  print "Reading big-endian binary files...\n";
	  $intform='N';
	  redo ENDIAN;
      }else{
	  die "Length of first record is $len instead of 8, 12, 16 or 24\n";
      }
  }
    open OUTFILE, ">$outdir/$rstfile";
    print OUTFILE $_;
    close OUTFILE;
}

open(INFILE, "$indir/restart.H") or die "cannot open $indir/restart.H\n";
$_=<INFILE>;
close INFILE;

open(OUTFILE, ">$outdir/restart.H") or die "cannot open $outdir/restart.H\n";
print OUTFILE $_;
close OUTFILE;

exit 0;

###############################################################################
sub convert4{
    for($i=0; $i<length(); $i+=4){
	substr($_,$i,4)=reverse substr($_,$i,4);
    }
}

###############################################################################
sub convert8{

    # initialize pointer for the string
    $i=0;

    while ( $i < length() ){
	# Get length of record
	$len = unpack($intform,substr($_,$i,4));
	
	# Check if length is reasonable
	die "At position $i record length $len is too large?!\n" 
	    if $i+$len > length();

	# Reverse leading 4 byte length marker
	$lenfixed = reverse(substr($_,$i,4));
	substr($_,$i,4)=$lenfixed;
	
	# Reverse 8 byte reals/integers
	for($j=$i+4; $j<$i+$len; $j+=8){
	    substr($_,$j,8)=reverse(substr($_,$j,8));
	}

	# Reverse trailing 4 byte length marker
	substr($_,$j,4)=$lenfixed;
	
	# Step to the next record
	$i = $j + 4;
    }
}


