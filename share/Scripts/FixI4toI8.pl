#!/usr/bin/perl
#^CFG COPYRIGHT UM

if($#ARGV != 0){
    print "
Purpose: add extra 0-s to 4 byte integers to make them 8 byte integers.

Typical usage: 

   FixI4toI8.pl octree.rst > octree_cray.rst

This should be done on a machine with the same endianness as the input
octree.rst file, but the script should not be run on the Cray itself,
because the Perl interpreter does not interpret long integers correctly 
on the Cray.

The endianness of the resulting octree_cray.rst file can be changed
with ConvertRestart if necessary.

";
    exit;
}

# No end of line record
undef $/;

# Define 4 byte zero for 8 byte integers in Cray
$zero=pack('L',0);

# Read the whole file into $_
$_=<>;

# Convert the file into a 4 byte integer array
@in=unpack 'L*', $_;
print STDERR "Number of 4 byte integers in file is ",$#in+1,"\n";

# Check if the length of the first record is reasonable
$first = $in[0];
die("Error: length of first record=$first is not 12\n") if $first!=12;

# Initialize array index
$i=0;

# Loop over the Fortran records
while($i <= $#in){

    # length of record in bytes
    $len=$in[$i];

    # number of 4 byte integers in record
    $n  = $len/4;

    # put 4 byte length at the beginning of record (twice the original!!!)
    push(@out,2*$len);

    # write 8 byte integers (0000xxxx) into @out
    for($j=$i+1;$j<=$i+$n;$j++){
	push(@out,$zero,$in[$j]);
    }

    # put 4 byte length at the end of record
    push(@out,2*$len);

    # next element
    $i=$i+$n+2;
}

$_=pack('L*',@out);

print;
