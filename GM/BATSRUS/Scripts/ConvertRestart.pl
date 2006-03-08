#!/usr/bin/perl
#^CFG COPYRIGHT UM

use strict;

#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Change Endianness of BATSRUS Restart Files with Scripts/ConvertRestart.pl}
#!ROUTINE: ConvertRestart.pl - change the endianness of BATSRUS restart files
#!DESCRIPTION:
# When moving restart files from one machine to another, it may be
# necessary to change the endianness of the binary files.
# This script is specifically written to convert all the binary restart
# files produced by BATSRUS, as well as to copy the ASCII header file.
#
#!REVISION HISTORY:
# 07/03/2001 G. Toth - initial revision
# 03/11/2003 G. Toth - generalized to work on a machine with 
#                      arbitrary endianness
# 03/08/2006 G. Toth - converts new restart file format and added use strict.
#EOP
&print_help if $#ARGV != 1;

my $indir =$ARGV[0];
my $outdir=$ARGV[1];

my $headerfile="restart.H";
my $datafile  ="data.rst";
my $octreefile="octree.rst";

my $ERROR = "ERROR in ConvertRestart.pl";

# Check things
opendir(INDIR, $indir) or die "$ERROR: cannot open input directory $indir!\n";

-d $outdir or mkdir($outdir,0777) or 
    die "$ERROR: cannot find/create output directory $outdir!\n";

$indir ne $outdir or 
    die "$ERROR: input and output directories must be different!\n";

# No end of line in files
undef $/;

# Copy restart.H from indir to outdir and obtain real precision
my $file = "$indir/$headerfile";
open(INFILE, $file) or die "$ERROR: cannot open $file\n"; 
$_=<INFILE>; 
close INFILE;

# Obtain the real precision
/([48])\s*nByteReal/i or die "$ERROR: could not find nByteReal in $file\n";
my $nByteReal = $1;

$file="$outdir/$headerfile";
open(OUTFILE, ">$file") or die "$ERROR: cannot open $file\n"; 
print OUTFILE $_;
close OUTFILE;

# There is 1 octree.rst file and many blk files
my @rstfiles = grep /^(blk\d*|octree|data)\.rst$/, readdir INDIR;
print "Number of files in $indir is ",1+$#rstfiles,"\n";

# Default integer format is VAX (little-endian, e.g. Linux PC)
my $intform='V';

print "Converting indir=$indir to outdir=$outdir...\n";
my $rstfile;
foreach $rstfile (@rstfiles){

    my $file="$indir/$rstfile";
    open(INFILE, $file) or die "$ERROR: cannot open $file!\n";
    $_=<INFILE>;
    close INFILE;

    # The $datafile is a binary direct access file with fixed record length.
    # Since there are no record markers it is trivial to swap the byte order.
    if($rstfile eq $datafile){
	if($nByteReal == 4){
	    &convert4;
	}else{
	    &convert8;
	}
    }else{
      ENDIAN:{
	  my $len = unpack $intform,substr($_,0,4);

	  if($len == 8 or $len == 12){
	      &convert4;
	  }elsif($len == 16 or $len == 24){
	      &convert_dbl;
	  }elsif($intform eq 'V'){
	      # Try Network integer format (big-endian, e.g. SGI, IBM, Cray)
	      print "Reading big-endian binary files...\n";
	      $intform='N';
	      redo ENDIAN;
	  }else{
	      die "$ERROR: Length of first record is $len ".
		  "instead of 8, 12, 16 or 24 in $file\n";
	  }
      }
    }
    my $file = "$outdir/$rstfile";
    open OUTFILE, ">$file" or die "$ERROR: cannot open $file\n";
    print OUTFILE $_;
    close OUTFILE;
}

exit 0;

###############################################################################
sub convert4{
    my $i;
    for($i=0; $i<length(); $i+=4){
	substr($_,$i,4)=reverse substr($_,$i,4);
    }
}

###############################################################################
sub convert8{
    my $i;
    for($i=0; $i<length(); $i+=8){
	substr($_,$i,8)=reverse substr($_,$i,8);
    }
}

###############################################################################
sub convert_dbl{

    # initialize pointer for the string
    my $i=0;

    while ( $i < length() ){
	# Get length of record
	my $len = unpack($intform,substr($_,$i,4));
	
	# Check if length is reasonable
	die "At position $i record length $len is too large?!\n" 
	    if $i+$len > length();

	# Reverse leading 4 byte length marker
	my $lenfixed = reverse(substr($_,$i,4));
	substr($_,$i,4)=$lenfixed;
	
	# Reverse 8 byte reals/integers
	my $j;
	for($j=$i+4; $j<$i+$len; $j+=8){
	    substr($_,$j,8)=reverse(substr($_,$j,8));
	}

	# Reverse trailing 4 byte length marker
	substr($_,$j,4)=$lenfixed;
	
	# Step to the next record
	$i = $j + 4;
    }
}

##############################################################################
sub print_help{

    print 
#BOC
"Purpose: change the endianness of *.rst files and copy restart.H.

Usage:

ConvertRestart.pl INPUTDIR OUTPUTDIR

This perl script can be run on a machine with arbitrary endianness.
The files are read from the INPUTDIR (e.g. restartOUT) and the result is 
saved into the OUTPUTDIR (e.g. restart_export). If the OUTPUTDIR does not
exist, it is created.

If the files are moved to/from a machine that uses 8 byte reals 
(some older Crays) then the octree file needs to be fixed for the 
4 vs. 8 byte integers using share/Scripts/FixI4toI8.pl or FixI8toI4.pl.
The endianness can be changed with ConvertRestart.pl both after or before the 
8-byte-integer-fix."
#EOC
    ,"\n\n";
    exit;
}
